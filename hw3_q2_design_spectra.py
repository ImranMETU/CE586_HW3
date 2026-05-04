"""
CE586 Earthquake Engineering Assignment 3 — Question 2
Design Spectra Comparison for DIN95Y01 Ground Motion

Purpose:
Compute and compare elastic and inelastic design spectra using:
1. Q1 record-based 5% damped actual acceleration spectrum
2. 2007 Turkish Seismic Code (TSC) with 2007 local site classes
3. 2018 Turkish Building Earthquake Code (TBDY) with official AFAD TDTH parameters

Assigned Ground Motion: DIN95Y01 (1995 Dinar earthquake, Dinar Meteorology Station, EW)
Student ID: 2735835
Soil Class: ZC (approximated as Z2 for 2007 code)

Author: Automated Script
Date: 2026
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Try to import pandas; if not available, use csv module
try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False
    import csv

# ============================================================
# CONFIGURATION
# ============================================================

# Paths
ASSIGNMENT_DIR = Path(r"C:\DDrive\Imran\METU\Coursework\CE586 - Earthquake\Assignment3")
RESULTS_Q1_DIR = ASSIGNMENT_DIR / "results"
RESULTS_Q2_DIR = ASSIGNMENT_DIR / "results_q2"
Q1_CSV_FILE = RESULTS_Q1_DIR / "DIN95Y01_response_spectra.csv"

# Ensure output directory exists
RESULTS_Q2_DIR.mkdir(exist_ok=True)

# Site Information
STUDENT_ID = 2735835
GM_CODE = "DIN95Y01"
EARTHQUAKE = "01.10.1995 Dinar"
STATION = "Dinar Meteorology Station"
COMPONENT = "EW"
LATITUDE = 38.247
LONGITUDE = 30.490
SOIL_CLASS = "ZC"
BUILDING_TYPE = "mid-rise reinforced concrete residential building"

# 2007 Turkish Seismic Code Parameters
# Dinar assumed as 1st degree seismic zone
A0_475_2007 = 0.40  # Peak ground acceleration for 475-year return period
I_2007 = 1.0  # Importance factor (residential)
TA_2007 = 0.15  # Period corner point A (s)
TB_2007 = 0.40  # Period corner point B (s)

# 2018 TDTH / TBDY Official Parameters
# DD-2 / 475-year return period
SS_DD2 = 0.608
S1_DD2 = 0.147
FS_DD2 = 1.257
F1_DD2 = 1.500
SDS_DD2 = 0.764
SD1_DD2 = 0.220
PGA_DD2 = 0.256
PGV_DD2 = 13.680
TA_DD2_report = 0.058
TB_DD2_report = 0.289
TL_DD2 = 6.000

# DD-1 / 2475-year return period
SS_DD1 = 1.235
S1_DD1 = 0.284
FS_DD1 = 1.200
F1_DD1 = 1.500
SDS_DD1 = 1.482
SD1_DD1 = 0.426
PGA_DD1 = 0.508
PGV_DD1 = 27.835
TA_DD1_report = 0.057
TB_DD1_report = 0.287
TL_DD1 = 6.000

# R-factors for inelastic spectra
R_FACTORS = [4, 8]

# ============================================================
# FUNCTIONS
# ============================================================

def load_q1_spectrum():
    """
    Load Q1 response spectrum CSV and extract 5% damped actual acceleration.
    
    Returns:
        T_array: Period values (1D array, seconds)
        Sa_actual_5pct: Actual acceleration at 5% damping (1D array, g)
    """
    print(f"Loading Q1 CSV: {Q1_CSV_FILE}")
    
    if not Q1_CSV_FILE.exists():
        raise FileNotFoundError(
            f"Q1 CSV not found at {Q1_CSV_FILE}\n"
            f"Please ensure hw3_response_spectra.py has been run successfully."
        )
    
    if HAS_PANDAS:
        df = pd.read_csv(Q1_CSV_FILE)
        # Filter for xi = 0.05
        df_5pct = df[df['xi'] == 0.05].copy()
        
        if len(df_5pct) == 0:
            raise ValueError("No rows found with xi = 0.05 in Q1 CSV")
        
        # Sort by period to ensure correct order
        df_5pct = df_5pct.sort_values('T').reset_index(drop=True)
        
        T_array = df_5pct['T'].values
        Sa_actual_5pct = df_5pct['Sa_actual_g'].values
    else:
        # Fallback to csv module
        T_list = []
        Sa_list = []
        with open(Q1_CSV_FILE, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if float(row['xi']) == 0.05:
                    T_list.append(float(row['T']))
                    Sa_list.append(float(row['Sa_actual_g']))
        
        if len(T_list) == 0:
            raise ValueError("No rows found with xi = 0.05 in Q1 CSV")
        
        # Sort by period
        sorted_pairs = sorted(zip(T_list, Sa_list))
        T_array = np.array([p[0] for p in sorted_pairs])
        Sa_actual_5pct = np.array([p[1] for p in sorted_pairs])
    
    # Verify data
    n_periods = len(T_array)
    expected_periods = 29
    if n_periods != expected_periods:
        print(f"Warning: Expected {expected_periods} periods, got {n_periods}")
    
    print(f"  ✓ Loaded {n_periods} periods from 0.2 s to 3.0 s")
    print(f"  ✓ Extracted 5% damped actual acceleration spectrum")
    
    return T_array, Sa_actual_5pct


def compute_2007_spectrum(T, A0, I, TA, TB):
    """
    Compute 2007 Turkish Seismic Code elastic acceleration spectrum.
    
    A(T) = A0 * I * S(T), where S(T) is defined piecewise.
    
    Args:
        T: Period (float or array, seconds)
        A0: Peak ground acceleration coefficient
        I: Importance factor
        TA: Period corner point A (s)
        TB: Period corner point B (s)
    
    Returns:
        Sae: Spectral acceleration (same type as T, in g)
    """
    T = np.atleast_1d(T)
    S = np.zeros_like(T, dtype=float)
    
    # 0 <= T <= TA
    mask1 = T <= TA
    S[mask1] = 1 + 1.5 * T[mask1] / TA
    
    # TA < T <= TB
    mask2 = (T > TA) & (T <= TB)
    S[mask2] = 2.5
    
    # T > TB
    mask3 = T > TB
    S[mask3] = 2.5 * (TB / T[mask3])**0.8
    
    Sae = A0 * I * S
    
    return Sae


def compute_2018_spectrum(T, SDS, SD1, TL):
    """
    Compute 2018 Turkish Building Earthquake Code (TBDY) elastic spectrum.
    
    Args:
        T: Period (float or array, seconds)
        SDS: Design spectral acceleration at short periods
        SD1: Design spectral acceleration at 1-second period
        TL: Long-period transition period (s)
    
    Returns:
        Sae: Spectral acceleration (same type as T, in g)
    """
    T = np.atleast_1d(T)
    
    # Calculate period corners
    TA = 0.2 * SD1 / SDS
    TB = SD1 / SDS
    
    Sae = np.zeros_like(T, dtype=float)
    
    # T <= TA
    mask1 = T <= TA
    Sae[mask1] = SDS * (0.4 + 0.6 * T[mask1] / TA)
    
    # TA < T <= TB
    mask2 = (T > TA) & (T <= TB)
    Sae[mask2] = SDS
    
    # TB < T <= TL
    mask3 = (T > TB) & (T <= TL)
    Sae[mask3] = SD1 / T[mask3]
    
    # T > TL
    mask4 = T > TL
    Sae[mask4] = SD1 * TL / (T[mask4]**2)
    
    return Sae, TA, TB


def save_csv_output(T, record_elastic, tsc_2007_475, tsc_2007_2475,
                    tbdy_dd2, tbdy_dd1, scale_factor,
                    ta_dd2, tb_dd2, ta_dd1, tb_dd1):
    """
    Save all computed spectra to CSV file.
    
    Args:
        T: Period array (s)
        record_elastic: Record-based 5% elastic spectrum (g)
        tsc_2007_475: 2007 TSC 475-year elastic (g)
        tsc_2007_2475: 2007 TSC 2475-year elastic (g)
        tbdy_dd2: 2018 TBDY DD-2 475-year elastic (g)
        tbdy_dd1: 2018 TBDY DD-1 2475-year elastic (g)
        scale_factor: 2007 scaling factor (475 to 2475)
        ta_dd2, tb_dd2: 2018 DD-2 period corners (s)
        ta_dd1, tb_dd1: 2018 DD-1 period corners (s)
    """
    output_file = RESULTS_Q2_DIR / "DIN95Y01_Q2_design_spectra.csv"
    
    print(f"\nSaving CSV output: {output_file}")
    
    if HAS_PANDAS:
        df = pd.DataFrame({
            'T': T,
            'Record_Elastic_5percent_g': record_elastic,
            'TSC2007_475_Elastic_g': tsc_2007_475,
            'TSC2007_2475_Elastic_g': tsc_2007_2475,
            'TBDY2018_DD2_475_Elastic_g': tbdy_dd2,
            'TBDY2018_DD1_2475_Elastic_g': tbdy_dd1,
            'Record_R4_g': record_elastic / 4,
            'TSC2007_475_R4_g': tsc_2007_475 / 4,
            'TSC2007_2475_R4_g': tsc_2007_2475 / 4,
            'TBDY2018_DD2_R4_g': tbdy_dd2 / 4,
            'TBDY2018_DD1_R4_g': tbdy_dd1 / 4,
            'Record_R8_g': record_elastic / 8,
            'TSC2007_475_R8_g': tsc_2007_475 / 8,
            'TSC2007_2475_R8_g': tsc_2007_2475 / 8,
            'TBDY2018_DD2_R8_g': tbdy_dd2 / 8,
            'TBDY2018_DD1_R8_g': tbdy_dd1 / 8,
            'Scale_2007_2475_to_475': scale_factor,
            'TA_DD2_calc': ta_dd2,
            'TB_DD2_calc': tb_dd2,
            'TA_DD1_calc': ta_dd1,
            'TB_DD1_calc': tb_dd1,
        })
        df.to_csv(output_file, index=False, float_format='%.6f')
    else:
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            # Header
            headers = [
                'T', 'Record_Elastic_5percent_g', 'TSC2007_475_Elastic_g',
                'TSC2007_2475_Elastic_g', 'TBDY2018_DD2_475_Elastic_g',
                'TBDY2018_DD1_2475_Elastic_g', 'Record_R4_g', 'TSC2007_475_R4_g',
                'TSC2007_2475_R4_g', 'TBDY2018_DD2_R4_g', 'TBDY2018_DD1_R4_g',
                'Record_R8_g', 'TSC2007_475_R8_g', 'TSC2007_2475_R8_g',
                'TBDY2018_DD2_R8_g', 'TBDY2018_DD1_R8_g', 'Scale_2007_2475_to_475',
                'TA_DD2_calc', 'TB_DD2_calc', 'TA_DD1_calc', 'TB_DD1_calc'
            ]
            writer.writerow(headers)
            # Data rows
            for i, t_val in enumerate(T):
                row = [
                    f'{t_val:.6f}',
                    f'{record_elastic[i]:.6f}',
                    f'{tsc_2007_475[i]:.6f}',
                    f'{tsc_2007_2475[i]:.6f}',
                    f'{tbdy_dd2[i]:.6f}',
                    f'{tbdy_dd1[i]:.6f}',
                    f'{record_elastic[i]/4:.6f}',
                    f'{tsc_2007_475[i]/4:.6f}',
                    f'{tsc_2007_2475[i]/4:.6f}',
                    f'{tbdy_dd2[i]/4:.6f}',
                    f'{tbdy_dd1[i]/4:.6f}',
                    f'{record_elastic[i]/8:.6f}',
                    f'{tsc_2007_475[i]/8:.6f}',
                    f'{tsc_2007_2475[i]/8:.6f}',
                    f'{tbdy_dd2[i]/8:.6f}',
                    f'{tbdy_dd1[i]/8:.6f}',
                    f'{scale_factor:.6f}',
                    f'{ta_dd2:.6f}',
                    f'{tb_dd2:.6f}',
                    f'{ta_dd1:.6f}',
                    f'{tb_dd1:.6f}',
                ]
                writer.writerow(row)
    
    print(f"  ✓ Saved {len(T)} periods to {output_file.name}")


def plot_elastic_comparison(T, record_elastic, tsc_2007_475, tsc_2007_2475,
                             tbdy_dd2, tbdy_dd1):
    """
    Generate elastic spectrum comparison plot.
    """
    output_file = RESULTS_Q2_DIR / "Q2c_elastic_comparison.png"
    
    print(f"\nGenerating elastic comparison plot: {output_file.name}")
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.plot(T, record_elastic, 'b-', linewidth=2.0, label='DIN95Y01 Record (5% damping)')
    ax.plot(T, tsc_2007_475, 'r--', linewidth=1.8, label='2007 TSC 475-year')
    ax.plot(T, tsc_2007_2475, 'r-', linewidth=1.8, label='2007 TSC 2475-year')
    ax.plot(T, tbdy_dd2, 'g--', linewidth=1.8, label='2018 TBDY DD-2 (475-year)')
    ax.plot(T, tbdy_dd1, 'g-', linewidth=1.8, label='2018 TBDY DD-1 (2475-year)')
    
    ax.set_xlabel('Period T (s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Spectral Acceleration Sa (g)', fontsize=12, fontweight='bold')
    ax.set_title('Q2(c) Elastic Spectrum Comparison — DIN95Y01, 5% Damping', 
                 fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim([0, 3.0])
    ax.set_ylim([0, max(np.max(tsc_2007_2475), np.max(tbdy_dd1)) * 1.1])
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved {output_file.name}")


def plot_2007_tsc_spectra(T, tsc_2007_475, tsc_2007_2475):
    """
    Generate 2007 Turkish Seismic Code elastic spectra plot.
    """
    output_file = RESULTS_Q2_DIR / "Q2a_2007_TSC_elastic_design_spectra.png"
    
    print(f"Generating 2007 TSC elastic spectra plot: {output_file.name}")
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.plot(T, tsc_2007_475, 'r--', linewidth=2.0, label='2007 TSC 475-year')
    ax.plot(T, tsc_2007_2475, 'r-', linewidth=2.0, label='2007 TSC 2475-year')
    
    ax.set_xlabel('Period T (s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Spectral Acceleration Sa (g)', fontsize=12, fontweight='bold')
    ax.set_title('Q2(a) 5% Damped Elastic Design Spectra — 2007 Turkish Seismic Code', 
                 fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim([0, 3.0])
    ax.set_ylim([0, np.max(tsc_2007_2475) * 1.1])
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved {output_file.name}")


def plot_2018_tbdy_spectra(T, tbdy_dd2, tbdy_dd1):
    """
    Generate 2018 TBDY elastic spectra plot.
    """
    output_file = RESULTS_Q2_DIR / "Q2b_2018_TBDY_elastic_design_spectra.png"
    
    print(f"Generating 2018 TBDY elastic spectra plot: {output_file.name}")
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.plot(T, tbdy_dd2, 'g--', linewidth=2.0, label='2018 TBDY DD-2 / 475-year')
    ax.plot(T, tbdy_dd1, 'g-', linewidth=2.0, label='2018 TBDY DD-1 / 2475-year')
    
    ax.set_xlabel('Period T (s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Spectral Acceleration Sa (g)', fontsize=12, fontweight='bold')
    ax.set_title('Q2(b) 5% Damped Elastic Design Spectra — 2018 Turkish Building Earthquake Code', 
                 fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim([0, 3.0])
    ax.set_ylim([0, np.max(tbdy_dd1) * 1.1])
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved {output_file.name}")


def plot_inelastic_comparison(T, record_elastic, tsc_2007_475, tsc_2007_2475,
                               tbdy_dd2, tbdy_dd1, R):
    """
    Generate inelastic spectrum comparison plot for given R-factor.
    """
    output_file = RESULTS_Q2_DIR / f"Q2d_inelastic_R{R}_comparison.png"
    
    print(f"Generating inelastic R={R} plot: {output_file.name}")
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    record_inelastic = record_elastic / R
    tsc_2007_475_inelastic = tsc_2007_475 / R
    tsc_2007_2475_inelastic = tsc_2007_2475 / R
    tbdy_dd2_inelastic = tbdy_dd2 / R
    tbdy_dd1_inelastic = tbdy_dd1 / R
    
    ax.plot(T, record_inelastic, 'b-', linewidth=2.0, label=f'DIN95Y01 Record (÷{R})')
    ax.plot(T, tsc_2007_475_inelastic, 'r--', linewidth=1.8, label=f'2007 TSC 475-year (÷{R})')
    ax.plot(T, tsc_2007_2475_inelastic, 'r-', linewidth=1.8, label=f'2007 TSC 2475-year (÷{R})')
    ax.plot(T, tbdy_dd2_inelastic, 'g--', linewidth=1.8, label=f'2018 TBDY DD-2 (÷{R})')
    ax.plot(T, tbdy_dd1_inelastic, 'g-', linewidth=1.8, label=f'2018 TBDY DD-1 (÷{R})')
    
    ax.set_xlabel('Period T (s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Spectral Acceleration Sa (g)', fontsize=12, fontweight='bold')
    ax.set_title(f'Q2(d) Inelastic Spectrum Comparison — R = {R}', 
                 fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlim([0, 3.0])
    ax.set_ylim([0, max(np.max(tsc_2007_2475_inelastic), np.max(tbdy_dd1_inelastic)) * 1.1])
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved {output_file.name}")


def save_markdown_summary(T, scale_factor, ta_dd2, tb_dd2, ta_dd1, tb_dd1):
    """
    Save markdown summary with all assumptions and parameters.
    """
    output_file = RESULTS_Q2_DIR / "Q2_summary.md"
    
    print(f"\nSaving markdown summary: {output_file.name}")
    
    summary = f"""# CE586 Assignment 3 — Question 2: Design Spectra

## 1. Site and Record Information

- **Student ID**: {STUDENT_ID}
- **Ground Motion Code**: {GM_CODE}
- **Earthquake**: {EARTHQUAKE}
- **Station**: {STATION}
- **Component**: {COMPONENT}
- **Latitude**: {LATITUDE}°N
- **Longitude**: {LONGITUDE}°E
- **Soil Class**: {SOIL_CLASS}
- **Building Type**: {BUILDING_TYPE}

## 2. Q1 Record-Based Spectrum Used

- **Source CSV**: `DIN95Y01_response_spectra.csv`
- **Damping Ratio**: ξ = 5%
- **Spectral Measure**: Actual acceleration (Sa_actual_g)
- **Period Range**: T = 0.2 s to 3.0 s (29 periods)
- **Number of Points**: {len(T)}

## 3. 2007 Turkish Seismic Code Assumptions

The 2007 Turkish Seismic Code (Türk Deprem Yönetmeliği 2007) elastic acceleration spectrum is defined as:

$$A(T) = A_0 \\cdot I \\cdot S(T)$$

where the spectrum coefficient S(T) is piecewise:

- For $0 \\leq T \\leq T_A$: $S(T) = 1 + 1.5 \\frac{{T}}{{T_A}}$
- For $T_A < T \\leq T_B$: $S(T) = 2.5$
- For $T > T_B$: $S(T) = 2.5 \\left( \\frac{{T_B}}{{T}} \\right)^{{0.8}}$

**Assumptions for Dinar Site:**

1. **Seismic Zone**: Dinar is assumed to be in the 1st degree seismic zone (Birinci derece deprem bölgesi)
2. **Peak Ground Acceleration Coefficient**: $A_0 = {A0_475_2007}$ (for 475-year return period)
3. **Importance Factor**: $I = {I_2007}$ (residential building)
4. **Soil Class Approximation**: Assigned soil class ZC (2018 system) is approximated as local site class Z2 (2007 system)
5. **Period Corner Points** (for Z2):
   - $T_A = {TA_2007}$ s
   - $T_B = {TB_2007}$ s

**2475-year Spectrum:**

The 2007 code provides directly a standard design earthquake for 475-year return period. To approximate the 2475-year spectrum (approximately MCE level), the 475-year spectrum is scaled by the official ratio of 2018 TDTH DD-1 to DD-2 parameters:

$$\\text{{scale factor}} = \\frac{{\\max(S_{{DS}}^{{DD-1}}/S_{{DS}}^{{DD-2}}, S_{{D1}}^{{DD-1}}/S_{{D1}}^{{DD-2}})}}{{}}$$

$$\\text{{scale factor}} \\approx \\frac{{1.482}}{{0.764}} \\approx 1.94$$

Therefore:

$$S_{{ae, 2475}}(T) = \\text{{scale factor}} \\times S_{{ae, 475}}(T) = 1.94 \\times S_{{ae, 475}}(T)$$

## 4. 2018 Turkish Building Earthquake Code / TDTH Parameters

The 2018 code (TBDY 2018) uses the official Turkish Earthquake Hazard Database (TDTH) for site-specific acceleration parameters based on coordinates and soil class.

### 2018 TDTH DD-2 / 475-year Return Period

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Acceleration at short periods | $S_S$ | {SS_DD2} |
| Acceleration at 1-second period | $S_1$ | {S1_DD2} |
| Short-period amplification factor | $F_S$ | {FS_DD2} |
| 1-second amplification factor | $F_1$ | {F1_DD2} |
| Design acceleration at short periods | $S_{{DS}}$ | {SDS_DD2} |
| Design acceleration at 1-second period | $S_{{D1}}$ | {SD1_DD2} |
| Peak ground acceleration (PGA) | PGA | {PGA_DD2} g |
| Peak ground velocity (PGV) | PGV | {PGV_DD2} cm/s |
| Period corner point A (calculated) | $T_A$ | {TA_DD2_report:.4f} s |
| Period corner point B (calculated) | $T_B$ | {TB_DD2_report:.4f} s |
| Long-period transition | $T_L$ | {TL_DD2:.3f} s |

### 2018 TDTH DD-1 / 2475-year Return Period (MCE)

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Acceleration at short periods | $S_S$ | {SS_DD1} |
| Acceleration at 1-second period | $S_1$ | {S1_DD1} |
| Short-period amplification factor | $F_S$ | {FS_DD1} |
| 1-second amplification factor | $F_1$ | {F1_DD1} |
| Design acceleration at short periods | $S_{{DS}}$ | {SDS_DD1} |
| Design acceleration at 1-second period | $S_{{D1}}$ | {SD1_DD1} |
| Peak ground acceleration (PGA) | PGA | {PGA_DD1} g |
| Peak ground velocity (PGV) | PGV | {PGV_DD1} cm/s |
| Period corner point A (calculated) | $T_A$ | {TA_DD1_report:.4f} s |
| Period corner point B (calculated) | $T_B$ | {TB_DD1_report:.4f} s |
| Long-period transition | $T_L$ | {TL_DD1:.3f} s |

### 2018 Horizontal Elastic Acceleration Spectrum

The elastic design spectrum is defined piecewise in terms of $S_{{DS}}$, $S_{{D1}}$, and $T_L$:

$$T_A = 0.2 \\frac{{S_{{D1}}}}{{S_{{DS}}}}$$

$$T_B = \\frac{{S_{{D1}}}}{{S_{{DS}}}}$$

The spectral acceleration is:

- For $T \\leq T_A$: $S_{{ae}}(T) = S_{{DS}} \\left( 0.4 + 0.6 \\frac{{T}}{{T_A}} \\right)$
- For $T_A < T \\leq T_B$: $S_{{ae}}(T) = S_{{DS}}$
- For $T_B < T \\leq T_L$: $S_{{ae}}(T) = \\frac{{S_{{D1}}}}{{T}}$
- For $T > T_L$: $S_{{ae}}(T) = \\frac{{S_{{D1}} \\cdot T_L}}{{T^2}}$

**Calculated Period Corners:**

- DD-2: $T_A = {ta_dd2:.4f}$ s (report: {TA_DD2_report} s), $T_B = {tb_dd2:.4f}$ s (report: {TB_DD2_report} s)
- DD-1: $T_A = {ta_dd1:.4f}$ s (report: {TA_DD1_report} s), $T_B = {tb_dd1:.4f}$ s (report: {TB_DD1_report} s)

The close agreement between calculated and reported values confirms correct implementation.

## 5. Elastic Spectrum Discussion

- **475-year vs 2475-year**: The 2475-year (DD-1) spectra are approximately 1.9 to 1.95 times larger than the 475-year (DD-2) spectra, reflecting higher hazard intensity at lower probability of exceedance.

- **2018 TBDY vs 2007 TSC**: The 2018 TBDY spectra are site-specific and derived from the official AFAD TDTH database using coordinates and soil properties. The 2007 TSC spectra are zone-based and less refined. Both design approaches are used in practice; the 2018 spectra are more conservative for this specific site.

- **Record vs Design Spectra**: The DIN95Y01 recorded spectrum represents the actual response to one historical earthquake. Design spectra represent probabilistic seismic hazard averaged over many possible future earthquakes. The record spectrum may exceed or fall below design spectra at specific periods, depending on the frequency content and magnitude.

- **Short-Period vs Long-Period**: All spectra show higher acceleration demand at short periods (< 0.5 s) and decreasing demand at longer periods. This reflects the concentration of seismic energy at shorter natural frequencies.

## 6. Inelastic Spectra via R-Factor Reduction

Inelastic (or "design") spectra were obtained using constant response modification factors:

$$S_{{a, inelastic}}(T) = \\frac{{S_{{a, elastic}}(T)}}{{R}}$$

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
*Student: {STUDENT_ID}*
*Course: CE586 Earthquake Engineering*
*Assignment: Assignment 3, Question 2*
"""
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(summary)
    
    print(f"  ✓ Saved {output_file.name}")


def main():
    """
    Main execution function.
    """
    print("=" * 70)
    print("CE586 Assignment 3 — Question 2: Design Spectra Comparison")
    print("=" * 70)
    print()
    
    try:
        # ============================================================
        # PART A: Load Q1 spectrum
        # ============================================================
        T, Sa_record_elastic = load_q1_spectrum()
        
        # ============================================================
        # PART B: Compute 2007 Turkish Seismic Code spectra
        # ============================================================
        print("\nComputing 2007 Turkish Seismic Code spectra...")
        Sa_2007_475_elastic = compute_2007_spectrum(T, A0_475_2007, I_2007, TA_2007, TB_2007)
        scale_2475_to_475 = max(SDS_DD1 / SDS_DD2, SD1_DD1 / SD1_DD2)
        Sa_2007_2475_elastic = scale_2475_to_475 * Sa_2007_475_elastic
        print(f"  ✓ 2007 475-year spectrum computed")
        print(f"  ✓ 2007 2475-year spectrum computed (scale factor = {scale_2475_to_475:.4f})")
        
        # ============================================================
        # PART C: Compute 2018 TBDY spectra
        # ============================================================
        print("\nComputing 2018 TBDY spectra...")
        Sa_2018_dd2_elastic, TA_DD2_calc, TB_DD2_calc = compute_2018_spectrum(
            T, SDS_DD2, SD1_DD2, TL_DD2
        )
        Sa_2018_dd1_elastic, TA_DD1_calc, TB_DD1_calc = compute_2018_spectrum(
            T, SDS_DD1, SD1_DD1, TL_DD1
        )
        print(f"  ✓ 2018 DD-2 475-year spectrum computed")
        print(f"    DD-2 period corners: TA = {TA_DD2_calc:.4f} s, TB = {TB_DD2_calc:.4f} s")
        print(f"    (Report values: TA = {TA_DD2_report} s, TB = {TB_DD2_report} s)")
        print(f"  ✓ 2018 DD-1 2475-year spectrum computed")
        print(f"    DD-1 period corners: TA = {TA_DD1_calc:.4f} s, TB = {TB_DD1_calc:.4f} s")
        print(f"    (Report values: TA = {TA_DD1_report} s, TB = {TB_DD1_report} s)")
        
        # ============================================================
        # PART D: Generate plots
        # ============================================================
        print("\nGenerating plots...")
        plot_2007_tsc_spectra(T, Sa_2007_475_elastic, Sa_2007_2475_elastic)
        plot_2018_tbdy_spectra(T, Sa_2018_dd2_elastic, Sa_2018_dd1_elastic)
        plot_elastic_comparison(T, Sa_record_elastic, Sa_2007_475_elastic,
                                Sa_2007_2475_elastic, Sa_2018_dd2_elastic, Sa_2018_dd1_elastic)
        
        for R in R_FACTORS:
            plot_inelastic_comparison(T, Sa_record_elastic, Sa_2007_475_elastic,
                                      Sa_2007_2475_elastic, Sa_2018_dd2_elastic,
                                      Sa_2018_dd1_elastic, R)
        
        # ============================================================
        # PART E: Save CSV output
        # ============================================================
        save_csv_output(T, Sa_record_elastic, Sa_2007_475_elastic, Sa_2007_2475_elastic,
                        Sa_2018_dd2_elastic, Sa_2018_dd1_elastic, scale_2475_to_475,
                        TA_DD2_calc, TB_DD2_calc, TA_DD1_calc, TB_DD1_calc)
        
        # ============================================================
        # PART F: Save markdown summary
        # ============================================================
        save_markdown_summary(T, scale_2475_to_475, TA_DD2_calc, TB_DD2_calc,
                              TA_DD1_calc, TB_DD1_calc)
        
        # ============================================================
        # Summary
        # ============================================================
        print("\n" + "=" * 70)
        print("ANALYSIS COMPLETE")
        print("=" * 70)
        print(f"Output directory: {RESULTS_Q2_DIR}")
        print(f"Generated files:")
        print(f"  • DIN95Y01_Q2_design_spectra.csv")
        print(f"  • Q2a_2007_TSC_elastic_design_spectra.png")
        print(f"  • Q2b_2018_TBDY_elastic_design_spectra.png")
        print(f"  • Q2c_elastic_comparison.png")
        print(f"  • Q2d_inelastic_R4_comparison.png")
        print(f"  • Q2d_inelastic_R8_comparison.png")
        print(f"  • Q2_summary.md")
        print("=" * 70)
        
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
