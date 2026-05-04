# CE586 Assignment 3 — Question 2: Design Spectra

## 1. Site and Record Information

- **Student ID**: 2735835
- **Ground Motion Code**: DIN95Y01
- **Earthquake**: 01.10.1995 Dinar
- **Station**: Dinar Meteorology Station
- **Component**: EW
- **Latitude**: 38.247°N
- **Longitude**: 30.49°E
- **Soil Class**: ZC
- **Building Type**: mid-rise reinforced concrete residential building

## 2. Q1 Record-Based Spectrum Used

- **Source CSV**: `DIN95Y01_response_spectra.csv`
- **Damping Ratio**: ξ = 5%
- **Spectral Measure**: Actual acceleration (Sa_actual_g)
- **Period Range**: T = 0.2 s to 3.0 s (29 periods)
- **Number of Points**: 29

## 3. 2007 Turkish Seismic Code Assumptions

The 2007 Turkish Seismic Code (Türk Deprem Yönetmeliği 2007) elastic acceleration spectrum is defined as:

$$A(T) = A_0 \cdot I \cdot S(T)$$

where the spectrum coefficient S(T) is piecewise:

- For $0 \leq T \leq T_A$: $S(T) = 1 + 1.5 \frac{T}{T_A}$
- For $T_A < T \leq T_B$: $S(T) = 2.5$
- For $T > T_B$: $S(T) = 2.5 \left( \frac{T_B}{T} \right)^{0.8}$

**Assumptions for Dinar Site:**

1. **Seismic Zone**: Dinar is assumed to be in the 1st degree seismic zone (Birinci derece deprem bölgesi)
2. **Peak Ground Acceleration Coefficient**: $A_0 = 0.4$ (for 475-year return period)
3. **Importance Factor**: $I = 1.0$ (residential building)
4. **Soil Class Approximation**: Assigned soil class ZC (2018 system) is approximated as local site class Z2 (2007 system)
5. **Period Corner Points** (for Z2):
   - $T_A = 0.15$ s
   - $T_B = 0.4$ s

**2475-year Spectrum:**

The 2007 code provides directly a standard design earthquake for 475-year return period. To approximate the 2475-year spectrum (approximately MCE level), the 475-year spectrum is scaled by the official ratio of 2018 TDTH DD-1 to DD-2 parameters:

$$\text{scale factor} = \frac{\max(S_{DS}^{DD-1}/S_{DS}^{DD-2}, S_{D1}^{DD-1}/S_{D1}^{DD-2})}{}$$

$$\text{scale factor} \approx \frac{1.482}{0.764} \approx 1.94$$

Therefore:

$$S_{ae, 2475}(T) = \text{scale factor} \times S_{ae, 475}(T) = 1.94 \times S_{ae, 475}(T)$$

## 4. 2018 Turkish Building Earthquake Code / TDTH Parameters

The 2018 code (TBDY 2018) uses the official Turkish Earthquake Hazard Database (TDTH) for site-specific acceleration parameters based on coordinates and soil class.

### 2018 TDTH DD-2 / 475-year Return Period

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Acceleration at short periods | $S_S$ | 0.608 |
| Acceleration at 1-second period | $S_1$ | 0.147 |
| Short-period amplification factor | $F_S$ | 1.257 |
| 1-second amplification factor | $F_1$ | 1.5 |
| Design acceleration at short periods | $S_{DS}$ | 0.764 |
| Design acceleration at 1-second period | $S_{D1}$ | 0.22 |
| Peak ground acceleration (PGA) | PGA | 0.256 g |
| Peak ground velocity (PGV) | PGV | 13.68 cm/s |
| Period corner point A (calculated) | $T_A$ | 0.0580 s |
| Period corner point B (calculated) | $T_B$ | 0.2890 s |
| Long-period transition | $T_L$ | 6.000 s |

### 2018 TDTH DD-1 / 2475-year Return Period (MCE)

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Acceleration at short periods | $S_S$ | 1.235 |
| Acceleration at 1-second period | $S_1$ | 0.284 |
| Short-period amplification factor | $F_S$ | 1.2 |
| 1-second amplification factor | $F_1$ | 1.5 |
| Design acceleration at short periods | $S_{DS}$ | 1.482 |
| Design acceleration at 1-second period | $S_{D1}$ | 0.426 |
| Peak ground acceleration (PGA) | PGA | 0.508 g |
| Peak ground velocity (PGV) | PGV | 27.835 cm/s |
| Period corner point A (calculated) | $T_A$ | 0.0570 s |
| Period corner point B (calculated) | $T_B$ | 0.2870 s |
| Long-period transition | $T_L$ | 6.000 s |

### 2018 Horizontal Elastic Acceleration Spectrum

The elastic design spectrum is defined piecewise in terms of $S_{DS}$, $S_{D1}$, and $T_L$:

$$T_A = 0.2 \frac{S_{D1}}{S_{DS}}$$

$$T_B = \frac{S_{D1}}{S_{DS}}$$

The spectral acceleration is:

- For $T \leq T_A$: $S_{ae}(T) = S_{DS} \left( 0.4 + 0.6 \frac{T}{T_A} \right)$
- For $T_A < T \leq T_B$: $S_{ae}(T) = S_{DS}$
- For $T_B < T \leq T_L$: $S_{ae}(T) = \frac{S_{D1}}{T}$
- For $T > T_L$: $S_{ae}(T) = \frac{S_{D1} \cdot T_L}{T^2}$

**Calculated Period Corners:**

- DD-2: $T_A = 0.0576$ s (report: 0.058 s), $T_B = 0.2880$ s (report: 0.289 s)
- DD-1: $T_A = 0.0575$ s (report: 0.057 s), $T_B = 0.2874$ s (report: 0.287 s)

The close agreement between calculated and reported values confirms correct implementation.

## 5. Elastic Spectrum Discussion

- **475-year vs 2475-year**: The 2475-year (DD-1) spectra are approximately 1.9 to 1.95 times larger than the 475-year (DD-2) spectra, reflecting higher hazard intensity at lower probability of exceedance.

- **2018 TBDY vs 2007 TSC**: The 2018 TBDY spectra are site-specific and derived from the official AFAD TDTH database using coordinates and soil properties. The 2007 TSC spectra are zone-based and less refined. Both design approaches are used in practice; the 2018 spectra are more conservative for this specific site.

- **Record vs Design Spectra**: The DIN95Y01 recorded spectrum represents the actual response to one historical earthquake. Design spectra represent probabilistic seismic hazard averaged over many possible future earthquakes. The record spectrum may exceed or fall below design spectra at specific periods, depending on the frequency content and magnitude.

- **Short-Period vs Long-Period**: All spectra show higher acceleration demand at short periods (< 0.5 s) and decreasing demand at longer periods. This reflects the concentration of seismic energy at shorter natural frequencies.

## 6. Inelastic Spectra via R-Factor Reduction

Inelastic (or "design") spectra were obtained using constant response modification factors:

$$S_{a, inelastic}(T) = \frac{S_{a, elastic}(T)}{R}$$

where:

- **R = 4**: Moderate ductility demand; suitable for ordinary moment-resisting frames or shear walls
- **R = 8**: High ductility demand; suitable for special moment-resisting frames or dual systems

This approach assumes that the structure can develop sufficient inelastic deformation capacity to dissipate the additional energy. Larger R factors reduce the design force (and thus construction cost and initial stiffness) but increase reliance on ductility.

**Key observations:**

- All inelastic spectra decrease proportionally with R: spectra for R = 8 are exactly half those for R = 4.
- Short-period structures are more sensitive to the choice of R because their elastic demands are highest.
- Long-period structures benefit more from ductility reduction because they are already in the declining portion of the spectrum.

## 7. Generated Files

All computed spectra and comparison plots have been saved to the `results_q2/` directory:

1. **DIN95Y01_Q2_design_spectra.csv**: Complete dataset with all elastic and inelastic spectra
2. **Q2a_2007_TSC_elastic_design_spectra.png**: 2007 TSC elastic spectra (475 and 2475 year)
3. **Q2b_2018_TBDY_elastic_design_spectra.png**: 2018 TBDY elastic spectra (DD-2 and DD-1)
4. **Q2c_elastic_comparison.png**: Elastic spectrum comparison (record vs 2007 TSC vs 2018 TBDY)
5. **Q2d_inelastic_R4_comparison.png**: Inelastic spectra with R = 4
6. **Q2d_inelastic_R8_comparison.png**: Inelastic spectra with R = 8
7. **Q2_summary.md**: This file

---

*Generated by: hw3_q2_design_spectra.py*
*Student: 2735835*
*Course: CE586 Earthquake Engineering*
*Assignment: Assignment 3, Question 2*
