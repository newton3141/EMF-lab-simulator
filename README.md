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
\[
\text{EMF} = - N \frac{d\Phi}{dt}
\]
\[
\Phi = \int B_z \, dA \approx A B_z(z)
\Rightarrow \text{EMF} = - N A \frac{dB_z}{dz} \frac{dz}{dt} = -N A \frac{dB_z}{dz} v
\]

### (2) 자석 자기장 근사 (유한 원통자석 모델)
\[
B_z(z) = \frac{B_r}{2}\left(
\frac{z+\frac{L}{2}}{\sqrt{(z+\frac{L}{2})^2 + R^2}}
-\frac{z-\frac{L}{2}}{\sqrt{(z-\frac{L}{2})^2 + R^2}}
\right)
\]
\[
\frac{dB_z}{dz} = \frac{B_r R^2}{2}
\left(
\frac{1}{(z+\frac{L}{2})^2+R^2)^{3/2}} - 
\frac{1}{(z-\frac{L}{2})^2+R^2)^{3/2}}
\right)
\]

### (3) 운동 방정식
\[
m\frac{dv}{dt} = mg - F_d - F_\text{lin} - F_\text{mag}
\]
\[
F_d = \frac{1}{2}\rho C_d A v|v|,\quad F_\text{lin} = k_\text{lin} v
\]
\[
\Rightarrow \frac{dv}{dt} = g - \frac{1}{2m}\rho C_d A v|v| - \frac{k_\text{lin}}{m}v - \frac{F_\text{mag}}{m}
\]

### (4) 회로 방정식 (Faraday–Lenz–Ohm law)
\[
L\frac{di}{dt} + Ri = \text{EMF}
\]
이를 암시적 오일러법(Implicit Euler) 으로 적분:
\[
i_{k} = \frac{i_{k-1} + \frac{\Delta t}{L}\text{EMF}_{k-1}}{1 + \frac{\Delta t}{L}R_\text{tot}}
\]

### (5) 렌츠힘 (자기 항력)
에너지 보존식에서 유도:
\[
P_\text{mech} = F_\text{mag}v = P_\text{elec} = \text{EMF}\cdot i
\Rightarrow F_\text{mag} = \frac{\text{EMF}\cdot i}{v}
\]

### (6) 코일 인덕턴스 (Wheeler 근사)
\[
L(\text{μH}) \approx \frac{N^2 D^2}{18D + 40\ell}
\]
\[
\text{또는} \quad L = 10^{-7} \frac{r^2 N^2}{9r + 10\ell} \quad [H]
\]

---

## 3. 설치 및 실행

### CLI 버전
```bash
pip install streamlit numpy pandas matplotlib
pip install streamlit numpy pandas matplotlib
streamlit run emf_lab_app_v2.py
