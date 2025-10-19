#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# EMF Lab Simulator (Pro) — v4 (foldered results + magnetic drag)

import argparse, math, os
from dataclasses import dataclass
from typing import Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RHO_CU = 1.68e-8

@dataclass
class Magnet:
    m: float
    R: float
    L: float
    Br: float

@dataclass
class Env:
    g: float = 9.80665
    rho: float = 1.2
    cd: float = 0.6
    k_lin: float = 0.0
    tube_R: float = 0.012

@dataclass
class Coil:
    N: int
    R: float
    L: float
    z: float
    name: str = ""
    wire_d: float = 0.0005
    R_override: float = None
    L_override: float = None
    def area(self): return math.pi*self.R**2
    def wire_length(self): return 2*math.pi*self.R*self.N
    def resistance(self):
        if self.R_override is not None: return self.R_override
        A = math.pi*(self.wire_d/2)**2
        return RHO_CU * self.wire_length() / max(A,1e-16)
    def inductance(self):
        if self.L_override is not None: return self.L_override
        r=self.R; l=max(self.L,1e-6)
        return 1e-7 * ((r**2)*(self.N**2)) / (9*r + 10*l)

def Bz_cyl(z_rel,R,L,Br):
    def u(x): return x/np.sqrt(x**2+R**2)
    return 0.5*Br*(u(z_rel+L/2)-u(z_rel-L/2))
def dBdz_cyl(z_rel,R,L,Br):
    def du(x): return (R**2)/np.power(x**2+R**2,1.5)
    return 0.5*Br*(du(z_rel+L/2)-du(z_rel-L/2))
def distributed(z_rel,R,L,Br,coilL,s=51):
    if coilL<=0: return Bz_cyl(z_rel,R,L,Br), dBdz_cyl(z_rel,R,L,Br)
    offs=np.linspace(-coilL/2,coilL/2,s)
    B=np.zeros_like(z_rel); dB=np.zeros_like(z_rel)
    for dz in offs:
        zr=z_rel+dz; B+=Bz_cyl(zr,R,L,Br); dB+=dBdz_cyl(zr,R,L,Br)
    return B/s, dB/s

def simulate(h0,z0,v0,tmax,dt,mag,env):
    n=int(np.ceil(tmax/dt))+1
    t=np.linspace(0,tmax,n); z=np.zeros_like(t); v=np.zeros_like(t)
    z[0]=z0; v[0]=v0; A=math.pi*env.tube_R**2
    for i in range(1,n):
        vi=v[i-1]; Fq=0.5*env.rho*env.cd*A*vi*abs(vi); Fl=env.k_lin*vi
        a=env.g-(Fq+Fl)/mag.m
        v[i]=vi+a*dt; z[i]=z[i-1]+v[i]*dt
        if z[i]>=h0: z[i:]=z[i]; v[i:]=0.0; break
    return t,z,v

def circuit(series: Dict[str,dict], mode="open", R_load=1e12, dt=5e-4, dt_sub=1):
    out={}; n=len(series["t"]); ds = dt/max(1,int(dt_sub))
    for label,d in series.items():
        if label in {"t","dt"}: continue
        e=d["emf"]; Lc=max(d["L"],1e-12); Rc=d["R"]
        if mode=="open":
            out[label]={"i":np.zeros_like(e),"v_meas":e}; continue
        Rtot = Rc + R_load
        i=np.zeros_like(e)
        for k in range(1,n):
            ik=i[k-1]; t_rem=dt
            while t_rem>1e-15:
                h = ds if t_rem>=ds else t_rem
                alpha = h/Lc
                ik = (ik + alpha*e[k-1]) / (1.0 + alpha*Rtot)
                t_rem -= h
            i[k]=ik
        out[label]={"i":i,"v_meas":R_load*i}
    return out

def magnetic_drag(series, circ, v):
    eps = 1e-12 + np.max(np.abs(v))*1e-12
    forces = {}
    n = len(v)
    F_total = np.zeros(n)
    for label in [k for k in circ.keys()]:
        e = series[label]["emf"]
        i = circ[label]["i"]
        F = np.zeros(n)
        nz = np.abs(v) > eps
        F[nz] = (e[nz] * i[nz]) / v[nz]
        forces[label] = F
        F_total += F
    return forces, F_total

def find_peaks(y,t,md=40):
    peaks=[]; last=-md
    for i in range(1,len(y)-1):
        if y[i]>y[i-1] and y[i]>y[i+1] and (i-last)>=md:
            peaks.append((t[i],y[i])); last=i
    return peaks

def make_save_dir(base="results"):
    import os
    os.makedirs(base, exist_ok=True)
    now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    path = os.path.join(base, now)
    orig = path; k=1
    while os.path.exists(path):
        path = f"{orig}_{k}"; k+=1
    os.makedirs(path, exist_ok=True)
    return path

def save_outputs(root,t,z,v,series,circ,F_forces,F_total):
    data={"t":t,"z":z,"v":v}
    for k,d in series.items():
        if k in {"t","dt"}: continue
        data[f"emf_{k}"]=d["emf"]; data[f"flux_{k}"]=d["flux"]; data[f"zrel_{k}"]=d["z_rel"]
    for k,d in circ.items():
        data[f"imeas_{k}"]=d["i"]; data[f"vmeas_{k}"]=d["v_meas"]
    for k,F in F_forces.items():
        data[f"Fmag_{k}"]=F
    data["Fmag_total"]=F_total
    pd.DataFrame(data).to_csv(root+"_timeseries.csv", index=False)

    rows=[]
    for k,d in series.items():
        if k in {"t","dt"}: continue
        for tp,yp in find_peaks(np.abs(d["emf"]), t, 40)[:2]:
            rows.append({"coil":k,"kind":"emf_abs_peak","t":tp,"val":yp})
    for k,d in circ.items():
        for tp,yp in find_peaks(np.abs(d["v_meas"]), t, 40)[:2]:
            rows.append({"coil":k,"kind":"v_meas_abs_peak","t":tp,"val":yp})
    pd.DataFrame(rows).to_csv(root+"_peaks.csv", index=False)

def plot_series(root,t,z,v,series,circ,coils,F_forces,F_total):
    plt.figure(figsize=(9,4.5)); plt.plot(t,z); plt.xlabel("time [s]"); plt.ylabel("z [m]"); plt.title("Magnet position"); plt.tight_layout(); plt.savefig(root+"_plot_z.png",dpi=150)
    plt.figure(figsize=(9,4.5)); plt.plot(t,v); plt.xlabel("time [s]"); plt.ylabel("v [m/s]"); plt.title("Magnet velocity"); plt.tight_layout(); plt.savefig(root+"_plot_v.png",dpi=150)
    plt.figure(figsize=(9,4.5))
    for k,d in series.items():
        if k in {"t","dt"}: continue
        plt.plot(t,d["emf"],label=k)
    plt.xlabel("time [s]"); plt.ylabel("emf [V] (model)"); plt.title("Induced EMF"); plt.legend(); plt.tight_layout(); plt.savefig(root+"_plot_emf.png",dpi=150)
    plt.figure(figsize=(9,4.5))
    for k,d in circ.items():
        plt.plot(t,d["v_meas"],label=k)
    plt.xlabel("time [s]"); plt.ylabel("V_load [V]"); plt.title("Measured voltage"); plt.legend(); plt.tight_layout(); plt.savefig(root+"_plot_vmeas.png",dpi=150)
    plt.figure(figsize=(9,4.5))
    for k,F in F_forces.items():
        plt.plot(t,F,label=f"Fmag_{k}")
    plt.plot(t,F_total,label="Fmag_total")
    plt.xlabel("time [s]"); plt.ylabel("F_mag [N]"); plt.title("Magnetic drag force (Lenz)"); plt.legend(); plt.tight_layout(); plt.savefig(root+"_plot_Fmag.png",dpi=150)
    labels=[k for k in series.keys() if k not in ("t","dt")]
    if labels:
        y0,y1=-0.015,0.015
        def z_to_x(zz): 
            return 0.5 if z.max()-z.min()<1e-12 else (zz - z.min())/(z.max()-z.min())
        first=labels[0]; import numpy as np
        idx=int(np.argmax(np.abs(series[first]["emf"])))
        plt.figure(figsize=(9,4.5))
        plt.plot([0,1],[y0,y0]); plt.plot([0,1],[y1,y1])
        for c in coils:
            xc=z_to_x(c.z); plt.plot([xc,xc],[y0*0.9,y1*0.9])
        xm=z_to_x(z[idx]); plt.scatter([xm],[0],s=200,marker="s")
        plt.xlim(-0.05,1.05); plt.ylim(y0*1.2,y1*1.2); plt.xlabel("tube (0 → bottom)"); plt.ylabel("radius")
        plt.title("Schematic at |emf| peak"); plt.tight_layout(); plt.savefig(root+"_plot_schematic.png",dpi=150)

def parse_coil(s: str) -> Coil:
    kw={}
    for part in s.split(","):
        if "=" in part:
            k,v=part.split("=",1); k=k.strip(); v=v.strip()
            if k in {"N"}: kw[k]=int(v)
            elif k in {"R","L","z","wire_d","wire"}:
                if k=="wire": k="wire_d"
                kw[k]=float(v)
            elif k=="name": kw[k]=v
            elif k=="R_override": kw[k]=float(v)
            elif k=="L_override": kw[k]=float(v)
    kw.setdefault("N",800); kw.setdefault("R",0.012); kw.setdefault("L",0.012); kw.setdefault("z",0.5)
    return Coil(**kw)

def run(params):
    t,z,v = simulate(params["h0"],params["z0"],params["v0"],params["t_max"],params["dt"],params["mag"],params["env"])
    series={"t":t,"dt":params["dt"]}
    for idx,c in enumerate(params["coils"]):
        zrel = z - c.z
        B,dB = distributed(zrel, params["mag"].R, params["mag"].L, params["mag"].Br, c.L)
        A = c.area(); emf = -c.N * A * dB * v
        label = c.name if c.name else f"coil_{idx+1}"
        series[label]={"emf":emf,"flux":c.N*A*B,"z_rel":zrel,"R":c.resistance(),"L":c.inductance()}
    circ = circuit(series, mode=params["mode"], R_load=params["R_load"], dt=params["dt"], dt_sub=params["dt_sub"])
    if params["mode"]!="open":
        F_forces, F_total = magnetic_drag(series, circ, v)
    else:
        zeros = np.zeros_like(t)
        F_forces = {k:np.zeros_like(series[k]["emf"]) for k in series if k not in {"t","dt"}}
        F_total = zeros
    return t,z,v,series,circ,F_forces,F_total

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--h0",type=float,default=0.9)
    ap.add_argument("--z0",type=float,default=0.0)
    ap.add_argument("--v0",type=float,default=0.0)
    ap.add_argument("--t_max",type=float,default=1.5)
    ap.add_argument("--dt",type=float,default=5e-4)
    ap.add_argument("--dt_sub",type=int,default=1, help="substeps per dt for circuit ODE (>=1)")
    ap.add_argument("--mag_m",type=float,default=0.025)
    ap.add_argument("--mag_R",type=float,default=0.006)
    ap.add_argument("--mag_L",type=float,default=0.010)
    ap.add_argument("--Br",type=float,default=1.2)
    ap.add_argument("--rho",type=float,default=1.2)
    ap.add_argument("--cd",type=float,default=0.6)
    ap.add_argument("--k_lin",type=float,default=0.0)
    ap.add_argument("--tube_R",type=float,default=0.012)
    ap.add_argument("--mode",choices=["open","series"],default="open")
    ap.add_argument("--R_load",type=float,default=1e12)
    ap.add_argument("--coil",action="append",default=[])
    ap.add_argument("--results_dir",type=str,default="results", help="root folder to save outputs")
    args=ap.parse_args()

    mag=Magnet(m=args.mag_m,R=args.mag_R,L=args.mag_L,Br=args.Br)
    env=Env(g=9.80665,rho=args.rho,cd=args.cd,k_lin=args.k_lin,tube_R=args.tube_R)
    coils=[parse_coil(s) for s in args.coil] if args.coil else [Coil(N=800,R=0.012,L=0.012,z=0.5,name="default")]
    params=dict(h0=args.h0,z0=args.z0,v0=args.v0,t_max=args.t_max,dt=args.dt,dt_sub=args.dt_sub,mag=mag,env=env,coils=coils,mode=args.mode,R_load=args.R_load)

    # results folder
    import os
    os.makedirs(args.results_dir, exist_ok=True)
    from datetime import datetime
    now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    save_dir = os.path.join(args.results_dir, now)
    orig = save_dir; k=1
    while os.path.exists(save_dir):
        save_dir = f"{orig}_{k}"; k+=1
    os.makedirs(save_dir, exist_ok=True)
    root = os.path.join(save_dir, "lab_run")

    t,z,v,series,circ,F_forces,F_total = run(params)
    save_outputs(root,t,z,v,series,circ,F_forces,F_total)
    plot_series(root,t,z,v,series,circ,coils,F_forces,F_total)

    print(f"Saved all outputs under: {save_dir}")

if __name__=="__main__":
    main()
