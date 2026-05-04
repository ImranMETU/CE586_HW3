# CE586 Assignment 3 Q1: Response Spectra Summary

## Assigned Ground Motion Record

| Parameter | Value |
|-----------|-------|
| **Student ID** | 2735835 |
| **GM Code** | DIN95Y01 |
| **Earthquake** | Dinar |
| **Date** | 01.10.1995 |
| **Station** | Dinar Meteorology Station |
| **Component** | EW |
| **Magnitude** | 6.0 |
| **Distance** | 5 km |
| **Latitude** | 38.247 |
| **Longitude** | 30.490 |
| **Soil Class** | ZC |

## Ground Motion Properties

| Property | Value |
|----------|-------|
| **Number of points** | 2796 |
| **Time step (dt)** | 0.0100 s |
| **Duration** | 27.96 s |
| **PGA (cm/s²)** | 313.190 |
| **PGA (g)** | 0.3194 |

## Analysis Parameters

### Newmark-Beta Method
- **Beta (β):** 0.25 (Constant Average Acceleration)
- **Gamma (γ):** 0.5
- **Unit mass (m):** 1.0 kg

### Response Spectra
- **Period range:** 0.2 – 3.0 s
- **Period increment:** 0.09999999999999998 s
- **Damping ratios (ξ):** 2%, 5%, 10%

## Formulations

### Equation of Motion
$$m \ddot{u} + c \dot{u} + k u = -m a_g(t)$$

where:
- $u$ = relative displacement
- $a_g(t)$ = ground acceleration
- $m$ = mass (1.0 kg)
- $c$ = damping coefficient = $2 \xi m \omega_n$
- $k$ = stiffness = $m \omega_n^2$
- $\omega_n$ = natural frequency = $2\pi / T$

### Displacement Spectrum
$$S_d = \max_t |u(t)|$$

### Pseudo-acceleration Spectrum
$$S_{a,pseudo} = \omega_n^2 S_d$$

### Actual Acceleration Spectrum
$$S_{a,actual} = \max_t |\ddot{u}(t) + a_g(t)|$$

## Generated Outputs

### Q1 Plots
1. **Q1a_displacement_spectra.png** – Displacement spectra for ξ = 2%, 5%, 10%
2. **Q1b_pseudo_acceleration_spectra.png** – Pseudo-acceleration spectra for ξ = 2%, 5%, 10%
3. **Q1c_actual_acceleration_spectra.png** – Actual acceleration spectra for ξ = 2%, 5%, 10%
4. **Q1d_pseudo_vs_actual_5percent.png** – Comparison of pseudo vs actual acceleration (ξ = 5%)

### Q2 Plots (Record-based Inelastic Spectra)
5. **Q2_record_inelastic_R4.png** – Inelastic spectrum with R = 4
6. **Q2_record_inelastic_R8.png** – Inelastic spectrum with R = 8

### Data Files
- **DIN95Y01_response_spectra.csv** – Complete tabulated results

## Discussion

### Effect of Damping on Response Spectra
Increasing damping reduces the spectral ordinates over most period ranges. This is 
expected because added damping dissipates energy, reducing the maximum response of the 
oscillator.

### Pseudo-acceleration vs Actual Acceleration
The pseudo-acceleration spectrum follows the general trend of the actual acceleration 
spectrum, but local differences occur because:

- **Pseudo-acceleration** is derived from the maximum relative displacement: 
  $S_{a,pseudo} = \omega_n^2 \max|u(t)|$
  
- **Actual spectral acceleration** is obtained from the absolute acceleration time history: 
  $S_{a,actual} = \max|\ddot{u}(t) + a_g(t)|$

The pseudo-acceleration assumes the maximum displacement and maximum acceleration occur 
at the same instant, which is generally not true. The actual acceleration spectrum captures 
the actual extremes of the total (absolute) acceleration response, making it more suitable 
for design of accelerating structures (e.g., equipment on flexible supports).

---
*Generated using Newmark-beta constant average acceleration method (β=1/4, γ=1/2)*
