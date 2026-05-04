"""
CE586 Earthquake Engineering - Assignment 3 Q1: Response Spectra Calculation
============================================================================

Computes elastic response spectra using the Newmark-beta method (constant average acceleration)
for the assigned ground motion record DIN95Y01 (Dinar, 1995-10-01).

Calculates:
  - Displacement spectra (Sd)
  - Pseudo-acceleration spectra (Sa,pseudo = omega_n² * Sd)
  - Actual absolute acceleration spectra (Sa,actual = max|a_rel + ag|)
  - Comparison of pseudo vs actual acceleration

Reference:
  - Equation of motion: m*u'' + c*u' + k*u = -m*ag(t)
  - Newmark-beta parameters: beta = 1/4, gamma = 1/2 (constant average acceleration)
  - Unit mass: m = 1.0
  
Student ID: 2735835
GM Code: DIN95Y01
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import warnings

warnings.filterwarnings("ignore", category=UserWarning)


# ============================================================================
# GLOBAL CONFIGURATION
# ============================================================================

BASE_DIR = Path(r"C:\DDrive\Imran\METU\Coursework\CE586 - Earthquake\Assignment3")
RESULTS_DIR = BASE_DIR / "results"
PARENT_DIR = BASE_DIR.parent

# Record metadata
RECORD_METADATA = {
    "Student ID": "2735835",
    "GM Code": "DIN95Y01",
    "Earthquake": "Dinar",
    "Date": "01.10.1995",
    "Station": "Dinar Meteorology Station",
    "Component": "EW",
    "Magnitude": "6.0",
    "Distance (km)": "5",
    "Latitude": "38.247",
    "Longitude": "30.490",
    "Soil Class": "ZC",
}

# Newmark-beta parameters
BETA = 1.0 / 4.0  # Constant average acceleration
GAMMA = 1.0 / 2.0
UNIT_MASS = 1.0

# Analysis parameters
PERIODS = np.round(np.arange(0.2, 3.0 + 0.1, 0.1), 1)
DAMPING_RATIOS = [0.02, 0.05, 0.10]


# ============================================================================
# NEWMARK-BETA SOLVER FOR GROUND MOTION EXCITATION
# ============================================================================

def newmark_ground_motion(ag, dt, T, xi, beta=0.25, gamma=0.5, m=1.0):
    """
    Solve the SDOF equation of motion using Newmark-beta method.
    
    Solves: m*u'' + c*u' + k*u = -m*ag(t)
    
    Parameters
    ----------
    ag : ndarray
        Ground acceleration time history (m/s²)
    dt : float
        Time step (seconds)
    T : float
        Natural period (seconds)
    xi : float
        Damping ratio (0.0 to 1.0 typically)
    beta : float
        Newmark-beta parameter (default 0.25 for constant average acceleration)
    gamma : float
        Newmark-gamma parameter (default 0.5 for constant average acceleration)
    m : float
        Mass (default 1.0)
    
    Returns
    -------
    u : ndarray
        Relative displacement (m)
    v : ndarray
        Relative velocity (m/s)
    a_rel : ndarray
        Relative acceleration (m/s²)
    """
    
    n_steps = len(ag)
    omega = 2.0 * np.pi / T
    k = m * omega**2
    c = 2.0 * xi * m * omega
    
    # Initialize response arrays
    u = np.zeros(n_steps)
    v = np.zeros(n_steps)
    a_rel = np.zeros(n_steps)
    
    # Initial conditions
    u[0] = 0.0
    v[0] = 0.0
    a_rel[0] = -ag[0]  # From m*a_rel[0] + 0 + 0 = -m*ag[0]
    
    # Newmark-beta coefficients (precomputed for efficiency)
    a0 = 1.0 / (beta * dt**2)
    a1 = gamma / (beta * dt)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2.0 * beta) - 1.0
    a4 = gamma / beta - 1.0
    a5 = dt * (gamma / (2.0 * beta) - 1.0)
    
    # Effective stiffness
    k_eff = k + a0 * m + a1 * c
    
    # Time-stepping loop
    for i in range(n_steps - 1):
        
        # Effective load: -m*ag(i+1) plus inertial and damping forces
        p_eff = (-m * ag[i + 1] + 
                 m * (a0 * u[i] + a2 * v[i] + a3 * a_rel[i]) +
                 c * (a1 * u[i] + a4 * v[i] + a5 * a_rel[i]))
        
        # Solve for displacement at next step
        u[i + 1] = p_eff / k_eff
        
        # Compute relative acceleration at next step
        a_rel[i + 1] = a0 * (u[i + 1] - u[i]) - a2 * v[i] - a3 * a_rel[i]
        
        # Compute relative velocity at next step
        v[i + 1] = v[i] + dt * ((1.0 - gamma) * a_rel[i] + gamma * a_rel[i + 1])
    
    return u, v, a_rel


# ============================================================================
# FILE HANDLING
# ============================================================================

def find_ground_motion_file():
    """
    Find the ground motion file automatically.
    
    Searches in Assignment3 directory and parent directory for files containing
    'DIN95Y01' with extensions .txt, .dat, .THF, or .csv.
    
    Returns
    -------
    Path
        Path to the ground motion file
    
    Raises
    ------
    FileNotFoundError
        If no suitable file is found
    """
    
    search_dirs = [BASE_DIR, PARENT_DIR]
    extensions = [".txt", ".dat", ".THF", ".csv"]
    
    for search_dir in search_dirs:
        if not search_dir.exists():
            continue
        
        for file in search_dir.rglob("*"):
            if "DIN95Y01" in file.name and file.suffix in extensions:
                print(f"✓ Found ground motion file: {file}")
                return file
    
    raise FileNotFoundError(
        f"Could not find DIN95Y01 ground motion file.\n"
        f"Searched in:\n"
        f"  - {BASE_DIR}\n"
        f"  - {PARENT_DIR}\n"
        f"Please place DIN95Y01.txt in {BASE_DIR}"
    )


def read_ground_motion_file(filepath):
    """
    Read ground motion file with two numeric columns.
    
    Parameters
    ----------
    filepath : Path
        Path to the ground motion file
    
    Returns
    -------
    time : ndarray
        Time array (seconds)
    ag_cm : ndarray
        Acceleration array (cm/s²)
    """
    
    data = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            try:
                parts = line.split()
                t = float(parts[0])
                a = float(parts[1])
                data.append([t, a])
            except (ValueError, IndexError):
                # Skip header or non-numeric lines
                continue
    
    if not data:
        raise ValueError(f"No numeric data found in {filepath}")
    
    data = np.array(data)
    time = data[:, 0]
    ag_cm = data[:, 1]
    
    return time, ag_cm


# ============================================================================
# RESPONSE SPECTRA COMPUTATION
# ============================================================================

def compute_response_spectra(time, ag_cm):
    """
    Compute elastic response spectra for multiple periods and damping ratios.
    
    Parameters
    ----------
    time : ndarray
        Time array (seconds)
    ag_cm : ndarray
        Ground acceleration array (cm/s²)
    
    Returns
    -------
    results_df : pd.DataFrame
        DataFrame with columns: T, xi, Sd_m, Sa_pseudo_mps2, Sa_pseudo_g, Sa_actual_mps2, Sa_actual_g
    """
    
    # Convert acceleration from cm/s² to m/s²
    ag = ag_cm / 100.0
    dt = time[1] - time[0]
    
    # Print ground motion information
    pga_cm = np.max(np.abs(ag_cm))
    pga_g = pga_cm / 980.665
    n_points = len(ag)
    
    print(f"\n{'─' * 70}")
    print("GROUND MOTION RECORD INFORMATION")
    print(f"{'─' * 70}")
    print(f"  Record: {RECORD_METADATA['GM Code']} ({RECORD_METADATA['Earthquake']})")
    print(f"  Date: {RECORD_METADATA['Date']}")
    print(f"  Station: {RECORD_METADATA['Station']}")
    print(f"  Component: {RECORD_METADATA['Component']}")
    print(f"  Number of points: {n_points}")
    print(f"  Time step (dt): {dt:.4f} s")
    print(f"  Duration: {time[-1]:.2f} s")
    print(f"  PGA: {pga_cm:.3f} cm/s² = {pga_g:.4f} g")
    print(f"{'─' * 70}\n")
    
    # Storage for results
    results = []
    
    print("Computing response spectra...")
    
    # Loop over damping ratios and periods
    for xi in DAMPING_RATIOS:
        print(f"  Damping ratio ξ = {xi*100:.1f}%", end=" ... ")
        
        for T in PERIODS:
            
            # Run Newmark solver
            u, v, a_rel = newmark_ground_motion(ag, dt, T, xi, 
                                                 beta=BETA, gamma=GAMMA, m=UNIT_MASS)
            
            # Compute spectral ordinates
            omega = 2.0 * np.pi / T
            
            # Displacement spectrum: max|u|
            Sd_m = np.max(np.abs(u))
            
            # Pseudo-acceleration spectrum: omega²*Sd
            Sa_pseudo_mps2 = omega**2 * Sd_m
            Sa_pseudo_g = Sa_pseudo_mps2 / 9.80665
            
            # Actual absolute acceleration: a_abs = a_rel + ag
            a_abs = a_rel + ag
            
            # Actual acceleration spectrum: max|a_abs|
            Sa_actual_mps2 = np.max(np.abs(a_abs))
            Sa_actual_g = Sa_actual_mps2 / 9.80665
            
            results.append({
                "T": T,
                "xi": xi,
                "Sd_m": Sd_m,
                "Sa_pseudo_mps2": Sa_pseudo_mps2,
                "Sa_pseudo_g": Sa_pseudo_g,
                "Sa_actual_mps2": Sa_actual_mps2,
                "Sa_actual_g": Sa_actual_g,
            })
        
        print("✓")
    
    results_df = pd.DataFrame(results)
    print(f"✓ Computed {len(results)} spectra\n")
    
    return results_df, ag, dt


def save_response_spectra_csv(results_df):
    """
    Save response spectra to CSV file.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Response spectra DataFrame
    """
    
    csv_path = RESULTS_DIR / "DIN95Y01_response_spectra.csv"
    results_df.to_csv(csv_path, index=False, float_format="%.6f")
    print(f"✓ Saved: {csv_path}")


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def plot_displacement_spectra(results_df):
    """
    Plot displacement spectra for all damping ratios.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Response spectra DataFrame
    """
    
    plt.figure(figsize=(12, 7))
    
    for xi in DAMPING_RATIOS:
        subset = results_df[results_df["xi"] == xi]
        plt.plot(subset["T"], subset["Sd_m"], marker="o", linewidth=2.5, 
                 markersize=6, label=f"ξ = {xi*100:.0f}%")
    
    plt.xlabel("Period T (s)", fontsize=13, fontweight="bold")
    plt.ylabel("Displacement Spectrum Sd (m)", fontsize=13, fontweight="bold")
    plt.title("Elastic Displacement Response Spectra\nDIN95Y01 (Dinar, EW Component)", 
              fontsize=14, fontweight="bold")
    plt.grid(True, alpha=0.3, linestyle="--")
    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    
    png_path = RESULTS_DIR / "Q1a_displacement_spectra.png"
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    print(f"✓ Saved: {png_path}")
    plt.close()


def plot_pseudo_acceleration_spectra(results_df):
    """
    Plot pseudo-acceleration spectra for all damping ratios.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Response spectra DataFrame
    """
    
    plt.figure(figsize=(12, 7))
    
    for xi in DAMPING_RATIOS:
        subset = results_df[results_df["xi"] == xi]
        plt.plot(subset["T"], subset["Sa_pseudo_g"], marker="o", linewidth=2.5,
                 markersize=6, label=f"ξ = {xi*100:.0f}%")
    
    plt.xlabel("Period T (s)", fontsize=13, fontweight="bold")
    plt.ylabel("Pseudo-acceleration Spectrum Sa,pseudo (g)", fontsize=13, fontweight="bold")
    plt.title("Elastic Pseudo-acceleration Response Spectra\nDIN95Y01 (Dinar, EW Component)",
              fontsize=14, fontweight="bold")
    plt.grid(True, alpha=0.3, linestyle="--")
    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    
    png_path = RESULTS_DIR / "Q1b_pseudo_acceleration_spectra.png"
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    print(f"✓ Saved: {png_path}")
    plt.close()


def plot_actual_acceleration_spectra(results_df):
    """
    Plot actual acceleration spectra for all damping ratios.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Response spectra DataFrame
    """
    
    plt.figure(figsize=(12, 7))
    
    for xi in DAMPING_RATIOS:
        subset = results_df[results_df["xi"] == xi]
        plt.plot(subset["T"], subset["Sa_actual_g"], marker="o", linewidth=2.5,
                 markersize=6, label=f"ξ = {xi*100:.0f}%")
    
    plt.xlabel("Period T (s)", fontsize=13, fontweight="bold")
    plt.ylabel("Actual Acceleration Spectrum Sa,actual (g)", fontsize=13, fontweight="bold")
    plt.title("Elastic Actual Acceleration Response Spectra\nDIN95Y01 (Dinar, EW Component)",
              fontsize=14, fontweight="bold")
    plt.grid(True, alpha=0.3, linestyle="--")
    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    
    png_path = RESULTS_DIR / "Q1c_actual_acceleration_spectra.png"
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    print(f"✓ Saved: {png_path}")
    plt.close()


def plot_pseudo_vs_actual(results_df):
    """
    Plot pseudo-acceleration vs actual acceleration for ξ = 5%.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Response spectra DataFrame
    """
    
    plt.figure(figsize=(12, 7))
    
    subset_5pct = results_df[results_df["xi"] == 0.05]
    
    plt.plot(subset_5pct["T"], subset_5pct["Sa_pseudo_g"], marker="o", linewidth=2.5,
             markersize=7, label="Pseudo-acceleration (ξ = 5%)", color="tab:blue")
    plt.plot(subset_5pct["T"], subset_5pct["Sa_actual_g"], marker="s", linewidth=2.5,
             markersize=7, label="Actual acceleration (ξ = 5%)", color="tab:red")
    
    plt.xlabel("Period T (s)", fontsize=13, fontweight="bold")
    plt.ylabel("Spectral Acceleration (g)", fontsize=13, fontweight="bold")
    plt.title("Pseudo-acceleration vs Actual Acceleration (ξ = 5%)\nDIN95Y01 (Dinar, EW Component)",
              fontsize=14, fontweight="bold")
    plt.grid(True, alpha=0.3, linestyle="--")
    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    
    png_path = RESULTS_DIR / "Q1d_pseudo_vs_actual_5percent.png"
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    print(f"✓ Saved: {png_path}")
    plt.close()


def plot_inelastic_spectra_q2(results_df, R_values=[4, 8]):
    """
    Plot inelastic spectra (record-based) for Q2.
    
    Uses the 5% actual acceleration spectrum as the elastic record spectrum.
    Inelastic spectrum is: Sa,inelastic = Sa,elastic / R
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Response spectra DataFrame
    R_values : list
        Ductility factors to plot
    """
    
    subset_5pct = results_df[results_df["xi"] == 0.05]
    
    # Plot for each ductility factor
    for R in R_values:
        plt.figure(figsize=(12, 7))
        
        # Elastic spectrum (ξ = 5%)
        plt.plot(subset_5pct["T"], subset_5pct["Sa_actual_g"], marker="o", linewidth=2.5,
                 markersize=7, label="Elastic (R = ∞, ξ = 5%)", color="tab:blue")
        
        # Inelastic spectrum
        Sa_inelastic = subset_5pct["Sa_actual_g"] / R
        plt.plot(subset_5pct["T"], Sa_inelastic, marker="s", linewidth=2.5,
                 markersize=7, label=f"Inelastic (R = {R}, ξ = 5%)", color="tab:red")
        
        plt.xlabel("Period T (s)", fontsize=13, fontweight="bold")
        plt.ylabel("Spectral Acceleration (g)", fontsize=13, fontweight="bold")
        plt.title(f"Record-based Inelastic Response Spectrum (R = {R})\nDIN95Y01 (Dinar, EW Component)",
                  fontsize=14, fontweight="bold")
        plt.grid(True, alpha=0.3, linestyle="--")
        plt.legend(fontsize=11, loc="best")
        plt.tight_layout()
        
        png_path = RESULTS_DIR / f"Q2_record_inelastic_R{R}.png"
        plt.savefig(png_path, dpi=150, bbox_inches="tight")
        print(f"✓ Saved: {png_path}")
        plt.close()


# ============================================================================
# SUMMARY REPORT
# ============================================================================

def create_summary_markdown(results_df, ag, dt):
    """
    Create a Markdown summary report.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Response spectra DataFrame
    ag : ndarray
        Ground acceleration (m/s²)
    dt : float
        Time step (seconds)
    """
    
    summary_text = f"""# CE586 Assignment 3 Q1: Response Spectra Summary

## Assigned Ground Motion Record

| Parameter | Value |
|-----------|-------|
| **Student ID** | {RECORD_METADATA['Student ID']} |
| **GM Code** | {RECORD_METADATA['GM Code']} |
| **Earthquake** | {RECORD_METADATA['Earthquake']} |
| **Date** | {RECORD_METADATA['Date']} |
| **Station** | {RECORD_METADATA['Station']} |
| **Component** | {RECORD_METADATA['Component']} |
| **Magnitude** | {RECORD_METADATA['Magnitude']} |
| **Distance** | {RECORD_METADATA['Distance (km)']} km |
| **Latitude** | {RECORD_METADATA['Latitude']} |
| **Longitude** | {RECORD_METADATA['Longitude']} |
| **Soil Class** | {RECORD_METADATA['Soil Class']} |

## Ground Motion Properties

| Property | Value |
|----------|-------|
| **Number of points** | {len(ag)} |
| **Time step (dt)** | {dt:.4f} s |
| **Duration** | {len(ag) * dt:.2f} s |
| **PGA (cm/s²)** | {np.max(np.abs(ag)) * 100:.3f} |
| **PGA (g)** | {np.max(np.abs(ag)) / 9.80665:.4f} |

## Analysis Parameters

### Newmark-Beta Method
- **Beta (β):** {BETA} (Constant Average Acceleration)
- **Gamma (γ):** {GAMMA}
- **Unit mass (m):** {UNIT_MASS} kg

### Response Spectra
- **Period range:** {PERIODS[0]} – {PERIODS[-1]} s
- **Period increment:** {PERIODS[1] - PERIODS[0]} s
- **Damping ratios (ξ):** 2%, 5%, 10%

## Formulations

### Equation of Motion
$$m \\ddot{{u}} + c \\dot{{u}} + k u = -m a_g(t)$$

where:
- $u$ = relative displacement
- $a_g(t)$ = ground acceleration
- $m$ = mass (1.0 kg)
- $c$ = damping coefficient = $2 \\xi m \\omega_n$
- $k$ = stiffness = $m \\omega_n^2$
- $\\omega_n$ = natural frequency = $2\\pi / T$

### Displacement Spectrum
$$S_d = \\max_t |u(t)|$$

### Pseudo-acceleration Spectrum
$$S_{{a,pseudo}} = \\omega_n^2 S_d$$

### Actual Acceleration Spectrum
$$S_{{a,actual}} = \\max_t |\\ddot{{u}}(t) + a_g(t)|$$

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
  $S_{{a,pseudo}} = \\omega_n^2 \\max|u(t)|$
  
- **Actual spectral acceleration** is obtained from the absolute acceleration time history: 
  $S_{{a,actual}} = \\max|\\ddot{{u}}(t) + a_g(t)|$

The pseudo-acceleration assumes the maximum displacement and maximum acceleration occur 
at the same instant, which is generally not true. The actual acceleration spectrum captures 
the actual extremes of the total (absolute) acceleration response, making it more suitable 
for design of accelerating structures (e.g., equipment on flexible supports).

---
*Generated using Newmark-beta constant average acceleration method (β=1/4, γ=1/2)*
"""
    
    summary_path = RESULTS_DIR / "HW3_summary.md"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(summary_text)
    
    print(f"✓ Saved: {summary_path}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main analysis workflow."""
    
    print("\n" + "=" * 70)
    print("CE586 ASSIGNMENT 3 Q1: RESPONSE SPECTRA CALCULATION")
    print("=" * 70)
    print(f"Record: {RECORD_METADATA['GM Code']} | Student ID: {RECORD_METADATA['Student ID']}")
    print("=" * 70 + "\n")
    
    # Create results directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"✓ Results directory: {RESULTS_DIR}\n")
    
    # Step 1: Find and read ground motion file
    print("Step 1: Loading ground motion file...")
    gm_file = find_ground_motion_file()
    
    time, ag_cm = read_ground_motion_file(gm_file)
    print(f"✓ Read {len(ag_cm)} acceleration points\n")
    
    # Step 2: Compute response spectra
    print("Step 2: Computing response spectra...")
    results_df, ag, dt = compute_response_spectra(time, ag_cm)
    
    # Step 3: Save results
    print("Step 3: Saving results...")
    save_response_spectra_csv(results_df)
    print()
    
    # Step 4: Generate plots
    print("Step 4: Generating plots...")
    plot_displacement_spectra(results_df)
    plot_pseudo_acceleration_spectra(results_df)
    plot_actual_acceleration_spectra(results_df)
    plot_pseudo_vs_actual(results_df)
    plot_inelastic_spectra_q2(results_df, R_values=[4, 8])
    print()
    
    # Step 5: Create summary
    print("Step 5: Creating summary report...")
    create_summary_markdown(results_df, ag, dt)
    print()
    
    print("=" * 70)
    print("✓ ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nResults saved in: {RESULTS_DIR}")
    print("\nGenerated files:")
    print("  - DIN95Y01_response_spectra.csv")
    print("  - Q1a_displacement_spectra.png")
    print("  - Q1b_pseudo_acceleration_spectra.png")
    print("  - Q1c_actual_acceleration_spectra.png")
    print("  - Q1d_pseudo_vs_actual_5percent.png")
    print("  - Q2_record_inelastic_R4.png")
    print("  - Q2_record_inelastic_R8.png")
    print("  - HW3_summary.md")
    print()


if __name__ == "__main__":
    main()
