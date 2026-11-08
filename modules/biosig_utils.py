"""Biosignature detection and analysis functions.
Code developed by Andrew Klingberg, adapted by Jake Lowe.
"""

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Detection threshold functions

def min_snr_for_prob(grid, SNR_list, R_list, P_req):
    """Find minimum SNR and corresponding R that achieves probability threshold."""
    minSNR = np.nan
    bestR = np.nan
    bestS = np.inf

    for r in range(grid.shape[0]):
        idx = np.where(grid[r, :] >= P_req)[0]
        if idx.size > 0:
            s = SNR_list[idx[0]]
            if s < bestS:
                bestS = s
                minSNR = s
                bestR = R_list[r]

    if np.isinf(bestS):
        return np.nan, np.nan

    return minSNR, bestR


def min_r_for_prob(grid, SNR_list, R_list, P_req):
    """Find minimum R and corresponding SNR that achieves probability threshold."""
    minR = np.nan
    bestSNR = np.nan
    bestRval = np.inf

    for s in range(grid.shape[1]):
        idx = np.where(grid[:, s] >= P_req)[0]
        if idx.size > 0:
            r = R_list[idx[0]]
            if r < bestRval:
                bestRval = r
                minR = r
                bestSNR = SNR_list[s]

    if np.isinf(bestRval):
        return np.nan, np.nan

    return minR, bestSNR


# Spectral processing

def blur_with_ref(lam, spec, Rtar, ref_wl):
    """Gaussian spectral blurring with resolving power Rtar."""
    fwhm = ref_wl / Rtar
    sig = fwhm / 2.355
    dlam = np.median(np.diff(lam))
    n = int(np.ceil(4 * sig / dlam))

    if n < 1:
        return spec.copy()

    x = np.arange(-n, n + 1) * dlam
    k = np.exp(-0.5 * (x / sig) ** 2)
    k /= np.sum(k)

    pad = np.concatenate([np.full(n, spec[0]), spec, np.full(n, spec[-1])])
    y = np.convolve(pad, k, mode="same")
    out = y[n:-n]

    return np.maximum(out, 0)


# Noise models

def add_noise_at_ref(lam, spec, targetSNR, ref_wl):
    """Add Gaussian noise using the spectrum at ref_wl to set the noise level."""
    iref = np.argmin(np.abs(lam - ref_wl))
    sig = np.abs(spec[iref]) / max(targetSNR, np.finfo(float).eps)
    noisy = spec + sig * np.random.randn(*spec.shape)

    return np.maximum(noisy, 0)


def add_noise_for_two_bands(lam, spec, targetSNR, red_nm, nir_nm, hw):
    """Add Gaussian noise using the NIR band to set the noise level."""
    iR = np.abs(lam - red_nm) <= hw
    iN = np.abs(lam - nir_nm) <= hw

    Fr = np.mean(spec[iR]) if np.any(iR) else 0.0
    Fn = np.mean(spec[iN]) if np.any(iN) else 0.0

    sig = Fn / max(targetSNR, np.finfo(float).eps)
    noisy = spec + sig * np.random.randn(*spec.shape)

    return np.maximum(noisy, 0.0), sig


# Feature metrics

def band_depth_with_snr(lam, spec, band):
    """Calculate absorption band depth and SNR."""
    iB = (lam >= band["band"][0]) & (lam <= band["band"][1])
    iL = (lam >= band["L"][0]) & (lam <= band["L"][1])
    iR = (lam >= band["R"][0]) & (lam <= band["R"][1])

    lamC = np.concatenate([lam[iL], lam[iR]])
    yC = np.concatenate([spec[iL], spec[iR]])

    if len(lamC) < 2:
        return 0.0, 0.0

    p = np.polyfit(lamC, yC, 1)
    cont = np.polyval(p, band["center"])
    yin = np.mean(spec[iB])
    depth = 1 - yin / max(cont, np.finfo(float).eps)

    yhat = np.polyval(p, lamC)
    resid = yC - yhat

    sig = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    sig = max(sig, 1e-8)
    snr = depth / (sig / max(cont, np.finfo(float).eps))

    return depth, snr


def ndvi_with_snr(lam, spec, red_nm, nir_nm, hw, min_ndvi=0.20, sig_injected=None):
    """Calculate NDVI and signal-to-noise ratio."""
    iR = np.abs(lam - red_nm) <= hw
    iN = np.abs(lam - nir_nm) <= hw

    Fr = np.mean(spec[iR]) if np.any(iR) else 0.0
    Fn = np.mean(spec[iN]) if np.any(iN) else 0.0

    denom = Fn + Fr
    if denom <= 0:
        return 0.0, 0.0

    ndvi = (Fn - Fr) / denom
    if ndvi < min_ndvi:
        return ndvi, 0.0

    nR = max(np.sum(iR), 1)
    nN = max(np.sum(iN), 1)

    if sig_injected is not None and sig_injected > 0:
        sig_Fr = sig_injected / np.sqrt(nR)
        sig_Fn = sig_injected / np.sqrt(nN)
    else:
        mask = iR | iN
        band_spec = spec[mask]
        if len(band_spec) >= 4:
            diffs = np.diff(band_spec)
            mad_diff = np.median(np.abs(diffs - np.median(diffs)))
            sig_channel = (1.4826 * mad_diff) / np.sqrt(2.0)
        else:
            sig_channel = 1e-8

        sig_channel = max(sig_channel, 1e-8)
        sig_Fr = sig_channel / np.sqrt(nR)
        sig_Fn = sig_channel / np.sqrt(nN)

    sig_ndvi = (2.0 / (denom**2)) * np.sqrt((Fr * sig_Fn)**2 + (Fn * sig_Fr)**2)
    sN = ndvi / max(sig_ndvi, 1e-8)

    return ndvi, sN


# Detection grid

def run_detection_grid(df_master, biomes, bands, filter_mode="poi",
                       SNR_list=None, R_list=None, detK=3,
                       VRE_red=665, VRE_NIR=865, VRE_halfwidth=10,
                       verbose=True):
    """Run the SNR/resolving-power detection grid for each region of interest."""

    if SNR_list is None:
        SNR_list = np.arange(5, 31, 2)
    if R_list is None:
        R_list = np.arange(60, 201, 10)

    wavelengths = np.array(df_master.columns[2:], dtype=float)

    all_band_names = list(bands.keys()) + ["VRE"]
    detection = {}

    start_time = time.time()

    for b in biomes:
        name = b["name"]
        w, s, e, n = b[filter_mode]

        spatial_mask = (
            (df_master["lon"] >= w) & (df_master["lon"] <= e) &
            (df_master["lat"] >= s) & (df_master["lat"] <= n)
        )

        df = df_master[spatial_mask]

        if len(df) == 0:
            print(f"Warning: No overlapping data coordinates found for {name}")
            continue

        X = df.iloc[:, 2:].to_numpy(dtype=float)
        nS = X.shape[0]

        detection[name] = {
            band_name: np.zeros((len(R_list), len(SNR_list)))
            for band_name in all_band_names
        }

        if verbose:
            print(f"Processing {name}: {nS} measurements...")

        for r_idx, Rtar in enumerate(R_list):
            for s_idx, SNRtar in enumerate(SNR_list):
                counts = {band_name: 0 for band_name in all_band_names}

                for i in range(nS):
                    spec0 = X[i, :]

                    for band_name, info in bands.items():
                        spec_b = blur_with_ref(
                            wavelengths, spec0, Rtar, info["center"]
                        )
                        spec_n = add_noise_at_ref(
                            wavelengths, spec_b, SNRtar, info["center"]
                        )
                        _, sigma_val = band_depth_with_snr(
                            wavelengths, spec_n, info
                        )

                        if sigma_val >= detK:
                            counts[band_name] += 1

                    spec_b_vre = blur_with_ref(
                        wavelengths, spec0, Rtar, VRE_NIR
                    )
                    spec_n_vre, sig_injected = add_noise_for_two_bands(
                        wavelengths, spec_b_vre, SNRtar,
                        VRE_red, VRE_NIR, VRE_halfwidth
                    )

                    _, sN = ndvi_with_snr(
                        wavelengths, spec_n_vre,
                        VRE_red, VRE_NIR, VRE_halfwidth,
                        sig_injected=sig_injected
                    )

                    if sN >= detK:
                        counts["VRE"] += 1

                for band_name in all_band_names:
                    detection[name][band_name][r_idx, s_idx] = counts[band_name] / nS

    if verbose:
        print(f"\nGrid search completed in {time.time() - start_time:.1f} seconds.")

    return detection


def print_detection_sanity_check(detection, verbose=True):
    """Print minimum and maximum detection fractions for each biome and band."""
    if verbose:
        print("\nGrid Search Sanity Check (Min / Max Fractions)")

    for biome in detection:
        print(f"\n{biome}:")
        for band_name, matrix_data in detection[biome].items():
            print(f"  -> {band_name.ljust(18)} min/max: {matrix_data.min():.4f} / {matrix_data.max():.4f}")


# SNR50 analysis

def find_snr50(snr_values, detection_values, threshold=0.5):
    """Find SNR at which detection fraction crosses threshold using linear interpolation."""
    snr_values = np.asarray(snr_values)
    detection_values = np.asarray(detection_values)

    order = np.argsort(snr_values)
    snr_values = snr_values[order]
    detection_values = detection_values[order]

    if detection_values[0] >= threshold:
        return snr_values[0]

    if detection_values[-1] < threshold:
        return np.nan

    idx = np.where(detection_values >= threshold)[0][0]

    if detection_values[idx] == threshold:
        return snr_values[idx]

    x1 = snr_values[idx - 1]
    x2 = snr_values[idx]
    y1 = detection_values[idx - 1]
    y2 = detection_values[idx]

    return x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)


def calculate_snr50(detection, SNR_list, R_list, detection_threshold=0.5):
    """Calculate SNR50 as a function of resolving power and summary statistics."""
    snr50_by_R = {}
    snr50_summary = {}
    snr50_mean = {}
    snr50_std = {}

    for biome in detection:
        snr50_by_R[biome] = {}
        snr50_summary[biome] = {}
        snr50_mean[biome] = {}
        snr50_std[biome] = {}

        for band_name in detection[biome]:
            values = []

            for r_idx, Rtar in enumerate(R_list):
                detection_vs_snr = detection[biome][band_name][r_idx, :]
                values.append(find_snr50(
                    SNR_list, detection_vs_snr, threshold=detection_threshold
                ))

            values = np.array(values)
            snr50_by_R[biome][band_name] = values

            valid_values = values[np.isfinite(values)]

            if len(valid_values) == 0:
                snr50_summary[biome][band_name] = np.nan
                snr50_mean[biome][band_name] = np.nan
                snr50_std[biome][band_name] = np.nan
            else:
                snr50_summary[biome][band_name] = np.median(valid_values)
                snr50_mean[biome][band_name] = np.mean(valid_values)
                snr50_std[biome][band_name] = np.std(valid_values, ddof=1) if len(valid_values) > 1 else 0.0

    snr50_df = pd.DataFrame(snr50_summary).T
    snr50_mean_df = pd.DataFrame(snr50_mean).T
    snr50_std_df = pd.DataFrame(snr50_std).T

    return snr50_by_R, snr50_df, snr50_mean_df, snr50_std_df


def get_snr50_stats(snr50_by_R):
    """Return median, min, max, and standard deviation of SNR50 across R."""
    snr50_stats = {}

    for biome in snr50_by_R:
        snr50_stats[biome] = {}

        for band_name, values in snr50_by_R[biome].items():
            values = values[np.isfinite(values)]

            if len(values) == 0:
                snr50_stats[biome][band_name] = {
                    "median": np.nan, "min": np.nan, "max": np.nan, "std": np.nan
                }
            else:
                snr50_stats[biome][band_name] = {
                    "median": np.median(values),
                    "min": np.min(values),
                    "max": np.max(values),
                    "std": np.std(values)
                }

    return snr50_stats


def print_snr50_stats(snr50_by_R):
    """Print SNR50 statistics across resolving powers."""
    snr50_stats = get_snr50_stats(snr50_by_R)

    for biome in snr50_stats:
        print(f"\n{biome}")

        for band_name, stats in snr50_stats[biome].items():
            print(
                f"{band_name:18s} "
                f"SNR50 = {stats['median']:.2f} "
                f"(range {stats['min']:.2f}–{stats['max']:.2f}, "
                f"σ={stats['std']:.2f})"
            )


# Plotting

def plot_snr50_heatmap(snr50_mean_df, snr50_std_df, figsize=(12, 4)):
    """Plot mean ± standard deviation SNR50 heatmap."""
    all_band_names = snr50_mean_df.columns.tolist()

    annot_df = snr50_mean_df.T.copy().astype(object)

    for band_name in all_band_names:
        for biome in snr50_mean_df.index:
            mean_val = snr50_mean_df.loc[biome, band_name]
            std_val = snr50_std_df.loc[biome, band_name]

            if np.isfinite(mean_val):
                annot_df.loc[band_name, biome] = f"{mean_val:.1f} ± {std_val:.1f}"
            else:
                annot_df.loc[band_name, biome] = "N/A"

    plt.figure(figsize=figsize)

    ax = sns.heatmap(
        snr50_mean_df.T,
        cmap=plt.cm.get_cmap("YlGn_r"),
        annot=annot_df,
        fmt="",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": r"Mean $SNR_{50}$"}
    )

    ax.set_xlabel("Region of Interest / Biome")
    ax.set_ylabel("Biosignature")
    ax.set_title(
        r"Mean $SNR_{50}$ for 50% Biosignature Detection "
        r"($R$ = 60–200)",
        fontsize=15,
        fontweight="bold"
    )

    plt.tight_layout()
    plt.show()


def plot_detection_heatmaps(detection, SNR_list, R_list, bands, figsize_per_panel=4):
    """Plot detection-fraction heatmaps with 50% threshold contours for each biome."""
    ncols = 5
    num_bands = len(next(iter(detection.values())))
    nrows = int(np.ceil(num_bands / ncols))
    all_band_names = list(next(iter(detection.values())).keys())

    for biome in detection:
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(figsize_per_panel * ncols, figsize_per_panel * nrows),
            sharex=False, sharey=False
        )
        axes = np.atleast_1d(axes).flatten()

        for idx, band_name in enumerate(all_band_names):
            ax = axes[idx]
            grid_data = detection[biome][band_name]
            df_grid = pd.DataFrame(grid_data, index=R_list, columns=SNR_list)

            sns.heatmap(
                df_grid.iloc[::-1],
                ax=ax,
                cmap="viridis",
                vmin=0.0,
                vmax=1.0,
                cbar_kws={"label": "Detection Fraction"} if (idx % ncols == ncols - 1) else {}
            )

            X, Y = np.meshgrid(
                np.arange(len(SNR_list)) + 0.5,
                np.arange(len(R_list)) + 0.5
            )

            if grid_data.min() < 0.5 < grid_data.max():
                ax.contour(
                    X, Y, grid_data[::-1, :],
                    levels=[0.5],
                    colors="red",
                    linestyles="--",
                    linewidths=2
                )

            if band_name == "VRE":
                ax.set_title("VRE (NDVI Ratio)", fontsize=10, fontweight="bold")
            else:
                ax.set_title(
                    f"{band_name} ({bands[band_name]['center']}nm)",
                    fontsize=10,
                    fontweight="bold"
                )

            if idx >= (nrows - 1) * ncols:
                ax.set_xlabel("Target $SNR$")
            else:
                ax.set_xlabel("")

            if idx % ncols == 0:
                ax.set_ylabel("Resolving Power ($R$)")
            else:
                ax.set_ylabel("")

        for blank_idx in range(num_bands, len(axes)):
            axes[blank_idx].axis("off")

        plt.suptitle(
            f"Biosignature Detection (50% Threshold Contour) — Biome: {biome}",
            fontsize=16,
            fontweight="bold",
            y=0.99
        )

        plt.tight_layout()
        plt.show()
