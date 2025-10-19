#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# EMF Lab Simulator (GUI): Streamlit app with full widget control and Run button

import math, os, io, shutil, zipfile
from dataclasses import dataclass
from typing import Dict, List
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

RHO_CU = 1.68e-8

# ------------------ Core models (same physics as v4) ------------------

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
    os.makedirs(base, exist_ok=True)
    now = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    path = os.path.join(base, now)
    orig = path; k=1
    while os.path.exists(path):
        path = f"{orig}_{k}"; k+=1
    os.makedirs(path, exist_ok=True)
    return path

def run_once(h0,z0,v0,t_max,dt,dt_sub,mag,env,coils,mode,R_load,results_root="results"):
    # simulate
    t,z,v = simulate(h0,z0,v0,t_max,dt,mag,env)
    series={"t":t,"dt":dt}
    for idx,c in enumerate(coils):
        zrel = z - c.z
        B,dB = distributed(zrel, mag.R, mag.L, mag.Br, c.L)
        A = c.area(); emf = -c.N*A*dB*v
        label = c.name if c.name else f"coil_{idx+1}"
        series[label]={"emf":emf,"flux":c.N*A*B,"z_rel":zrel,"R":c.resistance(),"L":c.inductance()}
    circ = circuit(series, mode=mode, R_load=R_load, dt=dt, dt_sub=dt_sub)
    F_forces, F_total = magnetic_drag(series, circ, v) if mode!="open" else ({k:np.zeros_like(series[k]["emf"]) for k in series if k not in {"t","dt"}}, np.zeros_like(t))

    # save folder
    save_dir = make_save_dir(results_root)
    root = os.path.join(save_dir, "lab_run")

    # save csv
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

    # plots
    plt.figure(figsize=(9,4.5)); plt.plot(t,z); plt.xlabel("time [s]"); plt.ylabel("z [m]"); plt.title("Magnet position"); plt.tight_layout(); plt.savefig(root+"_plot_z.png",dpi=150); plt.close()
    plt.figure(figsize=(9,4.5)); plt.plot(t,v); plt.xlabel("time [s]"); plt.ylabel("v [m/s]"); plt.title("Magnet velocity"); plt.tight_layout(); plt.savefig(root+"_plot_v.png",dpi=150); plt.close()
    plt.figure(figsize=(9,4.5))
    for k,d in series.items():
        if k in {"t","dt"}: continue
        plt.plot(t,d["emf"],label=k)
    plt.xlabel("time [s]"); plt.ylabel("emf [V] (model)"); plt.title("Induced EMF"); plt.legend(); plt.tight_layout(); plt.savefig(root+"_plot_emf.png",dpi=150); plt.close()
    plt.figure(figsize=(9,4.5))
    for k,d in circ.items():
        plt.plot(t,d["v_meas"],label=k)
    plt.xlabel("time [s]"); plt.ylabel("V_load [V]"); plt.title("Measured voltage"); plt.legend(); plt.tight_layout(); plt.savefig(root+"_plot_vmeas.png",dpi=150); plt.close()
    plt.figure(figsize=(9,4.5))
    for k,F in F_forces.items():
        plt.plot(t,F,label=f"Fmag_{k}")
    plt.plot(t,F_total,label="Fmag_total")
    plt.xlabel("time [s]"); plt.ylabel("F_mag [N]"); plt.title("Magnetic drag force (Lenz)"); plt.legend(); plt.tight_layout(); plt.savefig(root+"_plot_Fmag.png",dpi=150); plt.close()

    # schematic
    labels=[k for k in series.keys() if k not in ("t","dt")]
    if labels:
        y0,y1=-0.015,0.015
        def z_to_x(zz): 
            return 0.5 if z.max()-z.min()<1e-12 else (zz - z.min())/(z.max()-z.min())
        first=labels[0]
        idx=int(np.argmax(np.abs(series[first]["emf"])))
        plt.figure(figsize=(9,4.5))
        plt.plot([0,1],[y0,y0]); plt.plot([0,1],[y1,y1])
        for c in coils:
            xc=z_to_x(c.z); plt.plot([xc,xc],[y0*0.9,y1*0.9])
        xm=z_to_x(z[idx]); plt.scatter([xm],[0],s=200,marker="s")
        plt.xlim(-0.05,1.05); plt.ylim(y0*1.2,y1*1.2); plt.xlabel("tube (0 → bottom)"); plt.ylabel("radius")
        plt.title("Schematic at |emf| peak"); plt.tight_layout(); plt.savefig(root+"_plot_schematic.png",dpi=150); plt.close()

    return save_dir, root+"_timeseries.csv", root+"_plot_emf.png", root+"_plot_Fmag.png"

# ------------------ Streamlit UI ------------------

st.set_page_config(page_title="EMF Lab Simulator (GUI)", layout="wide")

st.title("🧲 EMF Lab Simulator (GUI)")
st.caption("낙하 자석–코일 유도 실험을 GUI로 제어하고, 결과를 폴더에 자동 저장합니다.")

with st.sidebar:
    st.header("Simulation")
    h0   = st.number_input("h0 [m]", value=0.9, min_value=0.05, step=0.05, format="%.3f")
    z0   = st.number_input("z0 [m]", value=0.0, format="%.3f")
    v0   = st.number_input("v0 [m/s]", value=0.0, format="%.3f")
    tmax = st.number_input("t_max [s]", value=1.5, min_value=0.1, step=0.1, format="%.2f")
    dt   = st.number_input("dt [s]", value=5e-4, min_value=1e-5, step=1e-5, format="%.5f")
    dt_sub = st.number_input("dt_sub [substeps]", value=1, min_value=1, step=1)

    st.header("Magnet")
    mag_m = st.number_input("mass m [kg]", value=0.025, min_value=0.001, step=0.001, format="%.3f")
    mag_R = st.number_input("radius R [m]", value=0.006, min_value=0.001, step=0.001, format="%.3f")
    mag_L = st.number_input("length L [m]", value=0.010, min_value=0.001, step=0.001, format="%.3f")
    Br    = st.number_input("remanence Br [T]", value=1.2, min_value=0.1, step=0.1, format="%.2f")

    st.header("Environment")
    rho   = st.number_input("air density ρ [kg/m³]", value=1.2, min_value=0.1, step=0.1, format="%.2f")
    cd    = st.number_input("drag coefficient Cd", value=0.6, min_value=0.0, step=0.1, format="%.2f")
    k_lin = st.number_input("linear drag k_lin [N·s/m]", value=0.0, min_value=0.0, step=0.01, format="%.2f")
    tube_R= st.number_input("tube radius [m]", value=0.012, min_value=0.001, step=0.001, format="%.3f")

    st.header("Circuit")
    mode  = st.selectbox("mode", ["open","series"], index=0)
    R_load= st.number_input("R_load [Ω]", value=1e12 if mode=="open" else 50.0, min_value=0.0, step=1.0, format="%.3f")

st.subheader("Coils")
st.caption("코일을 표 형태로 편집하세요. L_override/R_override는 빈 칸이면 자동 계산합니다.")
# Default coil table
default_coils = pd.DataFrame([
    {"name":"C1","N":800,"R":0.012,"L":0.012,"z":0.50,"wire_d":0.0005,"L_override":None,"R_override":None},
])
coil_df = st.data_editor(default_coils, num_rows="dynamic")
st.write("코일 개수:", len(coil_df))

colL, colR = st.columns(2)
with colL:
    results_root = st.text_input("Results root folder", value="results")
with colR:
    preset = st.selectbox("Quick preset", ["None","Light load (R=100Ω)","Heavier coil (N=1200)"], index=0)
    if preset == "Light load (R=100Ω)":
        mode = "series"; R_load = 100.0
    elif preset == "Heavier coil (N=1200)":
        coil_df.loc[:, "N"] = 1200

run_clicked = st.button("▶ Run simulation", type="primary")

if run_clicked:
    # Build objects
    mag = Magnet(m=mag_m, R=mag_R, L=mag_L, Br=Br)
    env = Env(g=9.80665, rho=rho, cd=cd, k_lin=k_lin, tube_R=tube_R)
    coils = []
    for _,row in coil_df.iterrows():
        coils.append(Coil(
            N=int(row["N"]), R=float(row["R"]), L=float(row["L"]), z=float(row["z"]),
            name=str(row["name"]), wire_d=float(row["wire_d"]),
            L_override=None if pd.isna(row["L_override"]) else float(row["L_override"]),
            R_override=None if pd.isna(row["R_override"]) else float(row["R_override"]),
        ))
    save_dir, csv_path, emf_png, fmag_png = run_once(h0,z0,v0,tmax,dt,dt_sub,mag,env,coils,mode,R_load,results_root)

    st.success(f"Saved all outputs under: {save_dir}")
    st.write("Timeseries CSV:", csv_path)
    st.image(emf_png, caption="Induced EMF")
    st.image(fmag_png, caption="Magnetic drag force")

    # offer ZIP download
    mem = io.BytesIO()
    with zipfile.ZipFile(mem, "w", zipfile.ZIP_DEFLATED) as zf:
        for fn in os.listdir(save_dir):
            full = os.path.join(save_dir, fn)
            zf.write(full, arcname=fn)
    mem.seek(0)
    st.download_button("⬇ Download run folder as ZIP", mem, file_name=os.path.basename(save_dir)+".zip")

st.markdown("---")
st.caption("Tip: mode='open'이면 자기 항력(Fmag)은 0으로 나옵니다. series 모드로 부하를 연결하면 Lenz drag가 계산됩니다.")