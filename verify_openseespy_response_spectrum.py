"""
OpenSeesPy Verification Script for CE586 Assignment 3 Q1
=========================================================

Computes elastic response spectrum using OpenSeesPy for damping ratio 5%.
Compares results against Newmark-beta implementation.

Record: DIN95Y01 (Dinar, EW component)
Periods: 0.2 to 3.0 s (0.1 s increment)
Damping: ξ = 5% (constant)

Ground motion: DIN95Y01.THF (time in s, acceleration in cm/s²)

This script uses OpenSeesPy's transient analysis with Newmark integration
(beta=0.5, gamma=0.25) and compares spectral ordinates against the reference
Newmark-beta implementation (beta=0.25, gamma=0.5).
"""

import numpy as np
import pandas as pd
from pathlib import Path
import sys
import warnings
import tempfile
import os

warnings.filterwarnings("ignore")

# ============================================================================
# OPENSEESPY IMPORT AND INSTALLATION HANDLING
# ============================================================================

OPENSEESPY_AVAILABLE = False
ops = None

try:
    import openseespy.opensees as ops
    OPENSEESPY_AVAILABLE = True
    print("[OK] OpenSeesPy is available")
except Exception as e:
    OPENSEESPY_AVAILABLE = False
    print(f"[!] OpenSeesPy not available: {str(e)[:60]}")


# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = Path(r"C:\DDrive\Imran\METU\Coursework\CE586 - Earthquake\Assignment3")
RESULTS_DIR = BASE_DIR / "results"
PARENT_DIR = BASE_DIR.parent

# Analysis parameters
XI_TARGET = 0.05  # 5% damping
PERIODS = np.round(np.arange(0.2, 3.0 + 0.1, 0.1), 1)
DT_ANALYSIS = 0.01  # seconds
G = 9.80665  # Gravitational acceleration


# ============================================================================
# FILE HANDLING
# ============================================================================

def find_and_read_ground_motion():
    """
    Find and read DIN95Y01.THF ground motion file.
    
    Returns
    -------
    time : ndarray
        Time array (seconds)
    ag_cm : ndarray
        Acceleration array (cm/s²)
    """
    
    search_dirs = [BASE_DIR, PARENT_DIR]
    extensions = [".thf", ".THF", ".txt", ".dat", ".csv"]
    
    gm_file = None
    for search_dir in search_dirs:
        if not search_dir.exists():
            continue
        for file in search_dir.rglob("*DIN95Y01*"):
            if file.suffix.lower() in [e.lower() for e in extensions]:
                gm_file = file
                break
        if gm_file:
            break
    
    if not gm_file:
        raise FileNotFoundError(
            f"Could not find DIN95Y01 ground motion file. "
            f"Searched in {BASE_DIR} and {PARENT_DIR}"
        )
    
    print(f"✓ Found: {gm_file}")
    
    # Read file
    data = []
    with open(gm_file, 'r') as f:
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
                continue
    
    if not data:
        raise ValueError(f"No numeric data found in {gm_file}")
    
    data = np.array(data)
    time = data[:, 0]
    ag_cm = data[:, 1]
    
    return time, ag_cm


# ============================================================================
# OPENSEESPY RESPONSE SPECTRUM COMPUTATION
# ============================================================================

def compute_openseespy_spectrum(time, ag_cm):
    """
    Compute response spectrum using OpenSeesPy with manual time stepping.
    
    Uses OpenSeesPy's Newmark integrator (gamma=0.5, beta=0.25)
    with uniform excitation pattern for ground motion.
    
    Parameters
    ----------
    time : ndarray
        Time array (seconds)
    ag_cm : ndarray
        Ground acceleration (cm/s²)
    
    Returns
    -------
    results_df : pd.DataFrame
        DataFrame with T, Sd_m, Sa_pseudo_g, Sa_actual_g
    """
    
    if not OPENSEESPY_AVAILABLE:
        print("✗ OpenSeesPy is not available. Cannot proceed.")
        return None
    
    # Convert acceleration to m/s²
    ag = ag_cm / 100.0
    dt_gm = time[1] - time[0]
    
    print(f"\n{'─' * 70}")
    print("OPENSEESPY VERIFICATION (Manual Time Stepping)")
    print(f"{'─' * 70}")
    print(f"Record: DIN95Y01")
    print(f"Time points: {len(ag)}")
    print(f"Duration: {time[-1]:.2f} s")
    print(f"Damping: ξ = {XI_TARGET*100:.1f}%")
    print(f"Integrator: Newmark (γ=0.5, β=0.25)\n")
    
    results = []
    
    for period_idx, T in enumerate(PERIODS):
        
        print(f"  T = {T:.1f} s ({period_idx+1}/{len(PERIODS)})", end=" ... ")
        
        try:
            # Reset model
            ops.wipeAll()
            ops.model('basic', '-ndm', 1, '-ndf', 1)
            
            # Nodes: 0 = base (fixed), 1 = mass
            ops.node(0, 0.0)
            ops.node(1, 0.0)
            ops.fix(0, 1)
            
            # Mass
            m = 1.0
            ops.mass(1, m)
            
            # Stiffness from period
            omega = 2.0 * np.pi / T
            k = m * omega**2
            
            # Spring element
            ops.uniaxialMaterial('Elastic', 1, k)
            ops.element('truss', 1, 0, 1, 1.0, 1)
            
            # Damping: Rayleigh with single mode
            alpha = 2.0 * XI_TARGET * omega
            ops.rayleigh(alpha, 0.0)
            
            # Analysis setup
            ops.constraints('Plain')
            ops.numberer('Plain')
            ops.system('BandGeneral')
            ops.test('NormDispIncr', 1.0e-8, 100, 0)
            ops.algorithm('Newton')
            ops.integrator('Newmark', 0.5, 0.25)  # constant average acceleration
            ops.analysis('Transient')
            
            # Time series for ground motion
            ops.timeSeries('Series', 1, '-time', *time, '-values', *ag)
            
            # Uniform excitation
            ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
            
            # Store responses
            u_history = []
            v_history = []
            a_history = []
            t_history = []
            
            # Manually step through analysis
            current_time = 0.0
            step = 0
            max_steps = int(time[-1] / DT_ANALYSIS) + 10
            
            while current_time < time[-1] and step < max_steps:
                # Run one analysis step
                ok = ops.analyze(1, DT_ANALYSIS)
                
                if ok != 0:
                    print(f"(step {step} failed)")
                    break
                
                # Get current time and responses
                current_time = ops.getTime()
                u = ops.nodeDisp(1, 1)  # node 1, DOF 1
                v = ops.nodeVel(1, 1)
                a = ops.nodeAccel(1, 1)
                
                u_history.append(u)
                v_history.append(v)
                a_history.append(a)
                t_history.append(current_time)
                
                step += 1
            
            # Convert to numpy
            u_history = np.array(u_history)
            v_history = np.array(v_history)
            a_history = np.array(a_history)
            t_history = np.array(t_history)
            
            if len(u_history) < 3:
                print("insufficient data")
                continue
            
            # Interpolate ground acceleration to analysis time steps
            ag_history = np.interp(t_history, time, ag)
            
            # Absolute acceleration
            a_abs_history = a_history + ag_history
            
            # Spectral values
            Sd_m = np.max(np.abs(u_history))
            Sa_pseudo_mps2 = omega**2 * Sd_m
            Sa_pseudo_g = Sa_pseudo_mps2 / G
            Sa_actual_mps2 = np.max(np.abs(a_abs_history))
            Sa_actual_g = Sa_actual_mps2 / G
            
            results.append({
                'T': T,
                'Sd_m': Sd_m,
                'Sa_pseudo_mps2': Sa_pseudo_mps2,
                'Sa_pseudo_g': Sa_pseudo_g,
                'Sa_actual_mps2': Sa_actual_mps2,
                'Sa_actual_g': Sa_actual_g,
            })
            
            print("✓")
            
        except Exception as e:
            print(f"✗ {str(e)[:40]}")
    
    results_df = pd.DataFrame(results)
    return results_df


# ============================================================================
# COMPARISON
# ============================================================================

def compare_results(newmark_csv_path):
    """
    Compare OpenSeesPy results against Newmark-beta implementation.
    
    Parameters
    ----------
    newmark_csv_path : Path
        Path to DIN95Y01_response_spectra.csv from Newmark implementation
    """
    
    print("\n" + "=" * 70)
    print("OPENSEESPY vs NEWMARK COMPARISON (ξ = 5%)")
    print("=" * 70)
    
    # Read Newmark results
    try:
        newmark_df = pd.read_csv(newmark_csv_path)
        newmark_5pct = newmark_df[newmark_df['xi'] == 0.05].copy()
    except FileNotFoundError:
        print(f"✗ Could not find: {newmark_csv_path}")
        return
    
    # Read OpenSeesPy results
    opensees_csv_path = RESULTS_DIR / "openseespy_verification_5percent.csv"
    try:
        opensees_df = pd.read_csv(opensees_csv_path)
    except FileNotFoundError:
        print(f"✗ Could not find: {opensees_csv_path}")
        return
    
    # Reset indices for alignment
    newmark_5pct = newmark_5pct.reset_index(drop=True)
    opensees_df = opensees_df.reset_index(drop=True)
    
    if len(newmark_5pct) != len(opensees_df):
        print(f"⚠ Row count mismatch: Newmark={len(newmark_5pct)}, OpenSees={len(opensees_df)}")
        return
    
    # Compare Sa_actual_g
    sa_actual_diff_abs = np.abs(newmark_5pct['Sa_actual_g'] - opensees_df['Sa_actual_g'])
    sa_actual_diff_pct = 100.0 * sa_actual_diff_abs / np.maximum(np.abs(newmark_5pct['Sa_actual_g']), 1e-6)
    
    # Compare Sd_m
    sd_diff_abs = np.abs(newmark_5pct['Sd_m'] - opensees_df['Sd_m'])
    sd_diff_pct = 100.0 * sd_diff_abs / np.maximum(np.abs(newmark_5pct['Sd_m']), 1e-6)
    
    print("\n" + "-" * 70)
    print("Sa_actual_g (Spectral Acceleration, g)")
    print("-" * 70)
    print(f"Max absolute difference: {np.max(sa_actual_diff_abs):.6f} g")
    print(f"Max relative difference: {np.max(sa_actual_diff_pct):.2f} %")
    print(f"Mean absolute difference: {np.mean(sa_actual_diff_abs):.6f} g")
    print(f"Mean relative difference: {np.mean(sa_actual_diff_pct):.2f} %")
    
    print("\n" + "-" * 70)
    print("Sd_m (Displacement Spectrum, m)")
    print("-" * 70)
    print(f"Max absolute difference: {np.max(sd_diff_abs):.6e} m")
    print(f"Max relative difference: {np.max(sd_diff_pct):.2f} %")
    print(f"Mean absolute difference: {np.mean(sd_diff_abs):.6e} m")
    print(f"Mean relative difference: {np.mean(sd_diff_pct):.2f} %")
    
    # Show period-by-period comparison
    print("\n" + "-" * 70)
    print("Period-by-Period Comparison")
    print("-" * 70)
    print(f"{'T (s)':>8} {'Newmark Sa':>14} {'OpenSees Sa':>14} {'Diff (%)':>12}")
    print("-" * 70)
    
    for idx in range(min(5, len(newmark_5pct))):  # Show first 5
        t_val = newmark_5pct.iloc[idx]['T']
        sa_nm = newmark_5pct.iloc[idx]['Sa_actual_g']
        sa_os = opensees_df.iloc[idx]['Sa_actual_g']
        pct_diff = 100.0 * abs(sa_nm - sa_os) / max(abs(sa_nm), 1e-6)
        print(f"{t_val:8.1f} {sa_nm:14.6f} {sa_os:14.6f} {pct_diff:12.2f}")
    
    print()


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main verification workflow."""
    
    print("\n" + "=" * 70)
    print("OPENSEESPY VERIFICATION SCRIPT")
    print("CE586 Assignment 3 Q1 - Response Spectra")
    print("=" * 70 + "\n")
    
    if not OPENSEESPY_AVAILABLE:
        print("✗ OpenSeesPy is required but not available.")
        print("\nTo install OpenSeesPy, run:")
        print("  pip install openseespy")
        print("\nOr in Conda:")
        print("  conda install -c conda-forge openseespy")
        return
    
    # Create results directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Read ground motion
    print("Step 1: Reading ground motion...")
    try:
        time, ag_cm = find_and_read_ground_motion()
        print(f"✓ Read {len(ag_cm)} points")
    except Exception as e:
        print(f"✗ Error: {e}")
        return
    
    # Step 2: Compute OpenSeesPy spectrum
    print("\nStep 2: Computing response spectrum...")
    opensees_df = compute_openseespy_spectrum(time, ag_cm)
    
    if opensees_df is None:
        print("✗ Failed to compute OpenSeesPy spectrum")
        return
    
    # Step 3: Save results
    print("\nStep 3: Saving results...")
    opensees_csv = RESULTS_DIR / "openseespy_verification_5percent.csv"
    opensees_df.to_csv(opensees_csv, index=False, float_format="%.6f")
    print(f"✓ Saved: {opensees_csv}")
    
    # Step 4: Compare
    print("\nStep 4: Comparing against Newmark results...")
    newmark_csv = RESULTS_DIR / "DIN95Y01_response_spectra.csv"
    compare_results(newmark_csv)
    
    print("\n" + "=" * 70)
    print("✓ VERIFICATION COMPLETE")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
