# EMF-lab-simulator

# EMF Lab Simulator (v4.1, 수행평가 제출용)
**Electronic Induction and Lenz Force Numerical Simulator**  
*Developer: Aiden Jung 정윤서 (Seoul Science High School, 2025)*

---

## 1. 개요

이 프로그램은 자석 낙하 실험(패러데이 전자기 유도 법칙 검증) 을 컴퓨터 상에서 수치적으로 재현하는 시뮬레이터이다.  
네오디뮴 자석이 솔레노이드 코일을 통과할 때 발생하는 유도기전력(EMF), 렌츠힘(Lenz Force), 속도 변화, 낙하 시간 등을 계산하여  
실제 실험에서 얻은 오실로스코프 파형과 비교·검증할 수 있다.

- 목적: 실험 오차의 원인을 정량적으로 분석  
- 활용: 수행평가 보고서, 오차 계산, 시각화 시뮬레이션  
- 특징
  - 현실 파라미터 기반 (코일, 자석, 공기저항, 부하저항 등)
  - 안정적인 암시적 회로 적분법(Implicit Euler)
  - 자동 폴더·시간명 저장, 렌츠힘 계산, mV 단위 지원
  - Streamlit GUI 인터페이스 지원

---

## 2. 물리 이론 요약

### (1) 패러데이의 전자기 유도 법칙
![수식](https://latex.codecogs.com/svg.latex?%5Ctext%7BEMF%7D%20%3D%20-N%20%5Cfrac%7Bd%5CPhi%7D%7Bdt%7D)
![수식](https://latex.codecogs.com/svg.latex?%5CPhi%20%3D%20%5Cint%20B_z%20%5C%2C%20dA%20%5Capprox%20A%20B_z%28z%29%20%5CRightarrow%20%5Ctext%7BEMF%7D%20%3D%20-N%20A%20%5Cfrac%7BdB_z%7D%7Bdz%7D%20v)

### (2) 자석 자기장 근사 (유한 원통자석 모델)
![수식](https://latex.codecogs.com/svg.latex?B_z%28z%29%20%3D%20%5Cfrac%7BB_r%7D%7B2%7D%5Cleft%28%5Cfrac%7Bz%2B%5Cfrac%7BL%7D%7B2%7D%7D%7B%5Csqrt%7B%28z%2B%5Cfrac%7BL%7D%7B2%7D%29%5E2%2BR%5E2%7D%7D%20-%20%5Cfrac%7Bz-%5Cfrac%7BL%7D%7B2%7D%7D%7B%5Csqrt%7B%28z-%5Cfrac%7BL%7D%7B2%7D%29%5E2%2BR%5E2%7D%7D%5Cright%29)

---
### (3) 운동 방정식
![수식](https://latex.codecogs.com/svg.latex?m%5Cfrac%7Bdv%7D%7Bdt%7D%20%3D%20mg%20-%20F_d%20-%20F_%7B%5Ctext%7Blin%7D%7D%20-%20F_%7B%5Ctext%7Bmag%7D%7D)
![수식](https://latex.codecogs.com/svg.latex?F_d%20%3D%20%5Cfrac%7B1%7D%7B2%7D%5Crho%20C_d%20A%20v%7Cv%7C%2C%20%5Cquad%20F_%7B%5Ctext%7Blin%7D%7D%20%3D%20k_%7B%5Ctext%7Blin%7D%7Dv)
### (4) 회로 방정식 (Faraday–Lenz–Ohm law)
![수식](https://latex.codecogs.com/svg.latex?L%5Cfrac%7Bdi%7D%7Bdt%7D%20%2B%20Ri%20%3D%20%5Ctext%7BEMF%7D)
![수식](https://latex.codecogs.com/svg.latex?i_%7Bk%7D%20%3D%20%5Cfrac%7Bi_%7Bk-1%7D%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7BL%7D%5Ctext%7BEMF%7D_%7Bk-1%7D%7D%7B1%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7BL%7DR_%7B%5Ctext%7Btot%7D%7D%7D)
### (5) 렌츠힘 (자기 항력)
에너지 보존식에서 유도:
![수식](https://latex.codecogs.com/svg.latex?P_%7B%5Ctext%7Bmech%7D%7D%20%3D%20F_%7B%5Ctext%7Bmag%7D%7Dv%20%3D%20P_%7B%5Ctext%7Belec%7D%7D%20%3D%20%5Ctext%7BEMF%7D%5Ccdot%20i)
![수식](https://latex.codecogs.com/svg.latex?F_%7B%5Ctext%7Bmag%7D%7D%20%3D%20%5Cfrac%7B%5Ctext%7BEMF%7D%5Ccdot%20i%7D%7Bv%7D)
### (6) 코일 인덕턴스 (Wheeler 근사)
![수식](https://latex.codecogs.com/svg.latex?L%28%5Cmu%20H%29%20%5Capprox%20%5Cfrac%7BN%5E2%20D%5E2%7D%7B18D%20%2B%2040%5Cell%7D)
---

## 3. 설치 및 실행

### CLI 버전
```bash
pip install streamlit numpy pandas matplotlib
pip install streamlit numpy pandas matplotlib
streamlit run emf_lab_app_v2.py
```

####인자(토글로 조작)

| **Category** | **Argument** | **Unit** | **Default** | **Description** |
|---------------|--------------|-----------|--------------|-----------------|
| **Simulation** | `--h0` | m | 0.9 | Drop height (tube length) |
|  | `--z0` | m | 0.0 | Initial position |
|  | `--v0` | m/s | 0.0 | Initial velocity |
|  | `--t_max` | s | 1.5 | Total simulation time |
|  | `--dt` | s | 5e-4 | Time step size |
|  | `--dt_sub` | - | 5 | Number of substeps for RL integration |
| **Magnet** | `--mag_m` | kg | 0.025 | Magnet mass |
|  | `--mag_R` | m | 0.006 | Magnet radius |
|  | `--mag_L` | m | 0.010 | Magnet length |
|  | `--Br` | T | 1.2 | Residual flux density (magnet strength) |
| **Environment** | `--rho` | kg/m³ | 1.2 | Air density |
|  | `--cd` | - | 0.6 | Drag coefficient |
|  | `--k_lin` | N·s/m | 0.0 | Linear viscous damping constant |
|  | `--tube_R` | m | 0.012 | Tube inner radius |
| **Circuit** | `--mode` | - | `series` | Circuit type (`open` or `series`) |
|  | `--R_load` | Ω | 50 | Load resistance |
|  | `--voltage_unit` | - | `mV` | Voltage output unit |
|  | `--results_dir` | - | `results` | Folder for results |
| **Coil** | `--coil` | - | - | Coil definition string (multiple allowed) |
---


출력되는 파일
---
lab_run_timeseries.csv		t, z, v, EMF, 전류, 전압, 자기항력 등 시계열 데이터
lab_run_peaks.csv		피크값 요약
lab_run_plot_z.png		자석 위치 vs 시간
lab_run_plot_v.png		자석 속도 vs 시간
lab_run_plot_emf.png		유도기전력 그래프
lab_run_plot_vmeas.png		부하전압 (V 또는 mV)
lab_run_plot_Fmag.png		자기항력(렌츠힘)
lab_run_plot_schematic.png	코일–자석 배치도


사용 자료 출처
---
진공 투자율	μ₀	4π×10⁻⁷ H/m	NIST CODATA (2019)
구리 비저항	ρ₍Cu₎	1.68×10⁻⁸ Ω·m	NIST Material Data
공기 밀도	ρ	1.2 kg/m³	20 °C, 해면 기준
중력 가속도	g	9.80665 m/s²	국제 중력 기준값
항력계수	C_d	0.5 ~ 0.7 (원통형)	NASA Glenn Research Center
네오디뮴 잔류자속밀도	B_r	1.0 ~ 1.48 T	Stanford Magnets, Arnold Magnetics
온도계수(Br)	α_T	–0.08 ~ –0.12 %/°C	Stanford Magnets
솔레노이드 인덕턴스	L	Wheeler 식	H.A. Wheeler, Proc. IRE (1928)

