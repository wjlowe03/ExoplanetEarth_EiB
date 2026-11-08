import os
import numpy as np
import pandas as pd
import xarray as xr
import earthaccess # type: ignore
import pickle
import matplotlib.pyplot as plt


def _process_radiance(path, nadir_idx=0):
    """Processes L1C file into a TOA reflectance dataset."""

    dt = xr.open_datatree(path)

    with xr.open_dataset(path) as root:
        d_earth_sun = float(root.attrs.get("sun_earth_distance", 1.0))

    obs_view = dt["observation_data"].ds.isel(number_of_views=nadir_idx)
    geo_view = dt["geolocation_data"].ds.isel(number_of_views=nadir_idx)
    view_angle_profile = dt["sensor_views_bands"].ds.isel(number_of_views=nadir_idx)
    toa_radiance = obs_view["i"]
    dim_axis_name = "intensity_bands_per_view"
    wv_var_name = "intensity_wavelength"
    wavelengths = view_angle_profile[wv_var_name]
    sza = geo_view["solar_zenith_angle"]

    radiance_ds = toa_radiance.to_dataset(name="radiance")
    sza_ds = sza.to_dataset(name="solar_zenith_angle")

    ds = xr.merge((radiance_ds, sza_ds, geo_view[["latitude", "longitude"]]))

    ds = ds.rename(
        {
            "latitude": "lat",
            "longitude": "lon",
            "solar_zenith_angle": "sza",
            "bins_along_track": "y",
            "bins_across_track": "x",
            dim_axis_name: "wavelength_3d",
        }
    )

    ds = ds.assign_coords(wavelength_3d=wavelengths.values)
    ds = ds.set_coords(["lat", "lon"])

    # Returns exactly ONE variable named 'data_var' so your downstream code works
    target = ds["radiance"].to_dataset(name="data_var")

    return target


_INSTRUMENT_REGISTRY = {
    "PACE_OCI_L1C_SCI": {
        "processor": _process_radiance,
        "granule_key": "granule",
        "mean_variable": "data_var",
    },
}


def extract_radiance_spectra(short_name, granule_dict, cache_dir=None, **kwargs):
    
    # Check registry that a processing helper is available
    if short_name not in _INSTRUMENT_REGISTRY:
        raise ValueError(f"Unsupported short name identifier: {short_name}")
    instrument = _INSTRUMENT_REGISTRY[short_name]
    item_name = granule_dict["name"]
    granule_search_str = granule_dict[instrument["granule_key"]]

    print(f"\n[Querying] Targeting granule match for '{item_name}'...")

    # Earthaccess search
    try:
        results = earthaccess.search_data(
            short_name=short_name,
            granule_name=granule_search_str,
        )

    except Exception:
        print("[NASA Server Error] CMR down or connection issue.")
        return None

    if not results:
        print(f"[Warning] No matching granules found for: {granule_search_str}")
        return None

    granule_name = results[0]["umm"]["GranuleUR"]

    print(f"[Match] Found: {granule_name}. Streaming content arrays...")

    # paths = earthaccess.open(results)

    # # Process according to instrument
    # ds = instrument["processor"](paths[0], **kwargs)
    # mean_var = instrument["mean_variable"]

    # # Spatial sampling
    # spectra_results = {
    #     "granule_mean": (ds[mean_var].mean(dim=["y", "x"], skipna=True)),
    #     "subsample_mean": None,
    #     "biome_mean": None,
    #     "ds_reference": ds,
    # }

    # print(f"Processing completed for '{item_name}'.")

    
    paths = earthaccess.open(results)

    try:
        # Process the remote granule.
        ds = instrument["processor"](paths[0], **kwargs)
        mean_var = instrument["mean_variable"]

        # CRITICAL: Materialize all data while the remote file is accessible.
        # This forces lazy xarray arrays to load into local memory.
        ds.load()

        # Detach the dataset from its underlying file handles.
        # The loaded arrays remain available in memory.
        ds.close()

        # Calculate summary statistics from the fully loaded dataset.
        spectra_results = {
            "granule_mean": ds[mean_var].mean(
                dim=["y", "x"], skipna=True
            ),
            "subsample_mean": None,
            "biome_mean": None,
            "ds_reference": ds,
        }

        print(f"Processing completed for '{item_name}'.")

        # Cache the fully materialized results.
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)

            safe_granule_id = (
                granule_name.split(":")[-1]
                if ":" in granule_name
                else granule_name
            )

            cache_filename = (
                f"{item_name}_{safe_granule_id}_cache.pkl"
            )
            cache_file = os.path.join(cache_dir, cache_filename)

            print(f"Saving processed results as {cache_filename}")

            with open(cache_file, "wb") as f:
                pickle.dump(
                    spectra_results,
                    f,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )

        return spectra_results

    finally:
        # Release the remote file handles regardless of success or failure.
        for path in paths:
            close = getattr(path, "close", None)
            if callable(close):
                close()

    # # Cache to a .pkl file
    # if cache_dir:
    #     os.makedirs(cache_dir, exist_ok=True)

    #     safe_granule_id = (
    #         granule_name.split(":")[-1]
    #         if ":" in granule_name
    #         else granule_name
    #     )

    #     cache_filename = (f"{item_name}_{safe_granule_id}_cache.pkl")
    #     cache_file = os.path.join(cache_dir, cache_filename)

    #     print(f"Saving processed results as {cache_filename}")

    #     with open(cache_file, "wb") as f:
    #         pickle.dump(spectra_results, f)

    # return spectra_results


def plot_radiance_comparison(data_dir):
    """
    Plotting:
    1. PACE OCI + Galileo NIMS Panel 1 comparison (~0.25-2.25 um; wavelength overlap)
    2. Galileo NIMS Panel 2 (~2.4-3.8 um)
    3. Galileo NIMS Panel 3 (~3.8-5.2 um)
    """

    # Load Galileo data
    files = ["rad1a.csv", "rad1b.csv", "rad1c.csv"]
    panels = [pd.read_csv(os.path.join(data_dir, fname)) for fname in files]

    # Convert Galileo radiance [(erg s^-1 cm^-2 sr^-1 um^-1) -> (W m^-2 sr^-1 um^-1)]:
    galileo_wavelengths = [panel["wl"].to_numpy(dtype=float) for panel in panels]
    galileo_radiances = [panel["rho"].to_numpy(dtype=float) * 1e-3 for panel in panels]

    # Load/filter PACE data:
    PACE_FILE = os.path.join(data_dir, "PA_Ocean2_PACE_OCI.20260613T233819.L1C.V3.5km.nc_cache.pkl")

    with open(PACE_FILE, "rb") as f:
        rf_data = pickle.load(f)

    rf_ds = rf_data["ds_reference"]
    rf_clean = rf_ds["data_var"].where(rf_ds["data_var"] != -32767) # Remove PACE fill value
    rf_w, rf_s, rf_e, rf_n = (-155.5, -6.0, -154.0, -5.0) # Pixels of interest
    rf_mask = ((rf_ds["lon"] >= rf_w) & (rf_ds["lon"] <= rf_e) & (rf_ds["lat"] >= rf_s) & (rf_ds["lat"] <= rf_n)) # Mask
    rf_poi_mean = rf_clean.where(rf_mask).mean(dim=["y", "x"], skipna=True)  # Spatial average

    # PACE wavelength (nm -> um):
    pacific_poi_wl = (rf_poi_mean["wavelength_3d"].values / 1000.0)
    pacific_poi_refl = rf_poi_mean.values

    # Plot
    galileo_color = "#5d6268"
    pace_color = "#4082d8"

    # Sagan et al. Fig. 1 Panel 1 w/ PACE comparison:
    fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=600)

    # PACE OCI
    ax.plot(
        pacific_poi_wl,
        pacific_poi_refl,
        color=pace_color,
        linestyle="-",
        linewidth=0.8,
        alpha=1.0,
        label="PACE OCI"
    )

    # Galileo NIMS Panel 1
    ax.plot(
        galileo_wavelengths[0],
        galileo_radiances[0],
        color=galileo_color,
        linestyle="-",
        linewidth=0.8,
        alpha=1.0,
        label="Galileo NIMS"
    )

    FONTSIZE = 8

    #ax.set_title("Earth Spectral Signatures: Modern PACE OCI vs. Galileo NIMS (1990)", fontsize=14, pad=12)
    #ax.set_xlabel("Wavelength (μm)", fontsize=FONTSIZE)
    #ax.set_ylabel(r"Spectral Radiance (W m$^{-2}$ sr$^{-1}$ μm$^{-1}$)", fontsize=FONTSIZE)
    #ax.legend(frameon=True, fontsize=6, loc="upper right")
    ax.tick_params(axis='both', labelsize=8)

    #NOTE: Uncomment to plot inset on Lowe et al. Fig. 1
    ax.set_xlim(0.68, 1.0) 
    ax.set_ylim(0, 10)

    plt.savefig("fig1_panel1.svg", format="svg", bbox_inches="tight")
    plt.tight_layout()
    plt.show()

    # Sagan et al. Fig. 1 Panels 2-3
    for i in [1, 2]:
        fig, ax = plt.subplots(figsize=(3.5, 1.75), dpi=600)

        FONTSIZE = 8

        ax.plot(
            galileo_wavelengths[i],
            galileo_radiances[i],
            color=galileo_color,
            linestyle="-",
            linewidth=0.8,
            alpha=1.0,
            label="Galileo NIMS"
        )

        ax.set_xlabel("Wavelength (μm)", fontsize=FONTSIZE)
        #ax.set_ylabel(r"Radiance (W m$^{-2}$ sr$^{-1}$ μm$^{-1}$)", fontsize=FONTSIZE)
        ax.tick_params(axis='both', labelsize=8)

        #plt.savefig(f"fig1_panel{i+1}.svg", format="svg", bbox_inches="tight")
        plt.tight_layout()
        plt.show()


def calculate_vmr(data_dir, gases):
    """Calculate PACE-derived VMRs for Sagan et al. (1993) comparison."""

    # Load PACE data
    PACE_FILE = os.path.join(data_dir, "PA_Ocean2_PACE_OCI.20260613T233819.L1C.V3.5km.nc_cache.pkl")

    with open(PACE_FILE, "rb") as f:
        rf_data = pickle.load(f)

    rf_ds = rf_data["ds_reference"]
    rf_clean = rf_ds["data_var"].where(rf_ds["data_var"] != -32767) # Remove PACE fill value
    rf_w, rf_s, rf_e, rf_n = (-155.5, -6.0, -154.0, -5.0) # Pixels of interest
    rf_mask = ((rf_ds["lon"] >= rf_w) & (rf_ds["lon"] <= rf_e) & (rf_ds["lat"] >= rf_s) & (rf_ds["lat"] <= rf_n)) # Mask

    # Volume mixing ratio (VMR) calculations
    # Open ocean, similar to cloudless region near Borneo referenced in Sagan et al. 1993
    poi_dataset = rf_clean.where(rf_mask, drop=True)

    # PACE wavelength (nm -> um):
    pacific_poi_wl = poi_dataset["wavelength_3d"].values / 1000.0
    pacific_poi_refl = poi_dataset.values

    # Grab N. Pacific data
    cached_data = {
        "wavelengths": pacific_poi_wl,
        "radiance": pacific_poi_refl,
        "solar_zenith": 30.0
    }

    # Constants (assumed by Sagan et al.)
    LOSCHMIDT = 2.68678e19        # molecules/cm^3 at STP
    TOTAL_AIR_COLUMN = 2.15e25     # molecules/cm^2 total vertical dry air column

    # Recompute Airmass factor based on your cached SZA
    sza_val = float(cached_data["solar_zenith"])
    eta = (1.0 / np.cos(np.radians(sza_val))) + 1.0

    # Extraction loop
    dim_name = "wavelength_3d"

    results = []
    print("Processing...")

    for gas, cfg in gases.items():

        bc, wl_l, wl_r = cfg["center"], cfg["L"], cfg["R"]
        sigma = cfg["cross_section"]

        # Extract matrix slices
        rad_l = poi_dataset.sel({dim_name: wl_l * 1000.0}, method="nearest")
        rad_r = poi_dataset.sel({dim_name: wl_r * 1000.0}, method="nearest")
        I_obs = poi_dataset.sel({dim_name: bc * 1000.0}, method="nearest")

        # Identify exact wavelengths matched by xarray to calculate the physical slope
        wl_actual_l = float(rad_l["wavelength_3d"].values / 1000.0)
        wl_actual_r = float(rad_r["wavelength_3d"].values / 1000.0)

        # Compute linear continuum baselines pixel-by-pixel across the grid
        if wl_actual_r == wl_actual_l:
            slope = 0.0
        else:
            slope = (rad_r - rad_l) / (wl_actual_r - wl_actual_l)

        I_0 = rad_l + slope * (bc - wl_actual_l)

        # Optical depth calculation
        ratio = I_obs / I_0
        ratio_clipped = np.clip(ratio, 1e-5, 0.99999)
        tau_2d = -np.log(ratio_clipped)

        # Absolute vertical column conversions
        u_slant_2d = tau_2d / sigma
        u_vert_2d = u_slant_2d / eta
        abundance_amagat_2d = u_vert_2d / LOSCHMIDT
        vmr_derived_2d = u_vert_2d / TOTAL_AIR_COLUMN

        # Final calculations
        mean_tau = float(np.nanmean(tau_2d))
        mean_amagat = float(np.nanmean(abundance_amagat_2d))
        mean_vmr = float(np.nanmean(vmr_derived_2d))
        spatial_std_vmr = float(np.nanstd(vmr_derived_2d))

        print(f"> Processed {gas:<4} | Mean VMR: {mean_vmr:.3e} | Spatial Std Dev: {spatial_std_vmr:.3e}")

        results.append({
            "Gas": gas,
            "Mean Optical Depth": f"{mean_tau:.4f}",
            "Mean Abundance (cm-amagat)": f"{mean_amagat:.4f}",
            "PACE Derived Mean VMR": f"{mean_vmr:.3e}",
            "PACE Spatial Std Dev": f"{spatial_std_vmr:.3e}"
        })

    df_sagan_summary = pd.DataFrame(results).set_index("Gas")

    return df_sagan_summary


def plot_reflectance_comparison(data_dir, correction_factor=None):
    """Plot PACE OCI and Galileo SSI spectral reflectance."""

    # Plot PACE/Galileo together:
    Fig3 = pd.read_csv(os.path.join(data_dir, "galileo/refl3.csv"))
    AreaA = os.path.join(data_dir, "toa_cache/Atacama_PACE_OCI.20260511T173132.L1C.V3.5km.nc_cache.pkl")
    AreaB = os.path.join(data_dir, "toa_cache/F3_B_PACE_OCI.20260628T172425.L1C.V3.5km.nc_cache.pkl")
    AreaC = os.path.join(data_dir, "toa_cache/Amazon_PACE_OCI_L1C_SCI_PACE_OCI.20250827T170528.L1C.V3.5km.nc_3_cache.pkl")

    # Bounding boxes for pixels of interest
    at_w, at_s, at_e, at_n = -65.5, -21, -64.5, -20 #-69.5, -23.5, -69.0, -23.0
    ar_w, ar_s, ar_e, ar_n = -59.5, -37.5, -59.0 , -37.0
    rf_w, rf_s, rf_e, rf_n = -65.0, -10.0, -60.0, -5.0

    # ATACAMA (Area A)
    with open(AreaA, "rb") as f:
        at_data = pickle.load(f)
    at_ds = at_data["ds_reference"]
    at_clean = at_ds["data_var"].where(at_ds["data_var"] != -32767)
    at_mask = (at_ds["lon"] >= at_w) & (at_ds["lon"] <= at_e) & (at_ds["lat"] >= at_s) & (at_ds["lat"] <= at_n)
    at_poi_mean = at_clean.where(at_mask).mean(dim=["y", "x"], skipna=True)
    atacama_poi_wl = at_poi_mean["wavelength_3d"].values / 1000.0
    atacama_poi_refl = at_poi_mean.values

    # ARGENTINA (Area B)
    with open(AreaB, "rb") as f:
        ar_data = pickle.load(f)
    ar_ds = ar_data["ds_reference"]
    ar_clean = ar_ds["data_var"].where(ar_ds["data_var"] != -32767)
    ar_mask = (ar_ds["lon"] >= ar_w) & (ar_ds["lon"] <= ar_e) & (ar_ds["lat"] >= ar_s) & (ar_ds["lat"] <= ar_n)
    ar_poi_mean = ar_clean.where(ar_mask).mean(dim=["y", "x"], skipna=True)
    argentina_poi_wl = ar_poi_mean["wavelength_3d"].values / 1000.0
    argentina_poi_refl = ar_poi_mean.values

    # RAINFOREST (Area C)
    with open(AreaC, "rb") as f:
        rf_data = pickle.load(f)
    rf_ds = rf_data["ds_reference"]
    rf_clean = rf_ds["data_var"].where(rf_ds["data_var"] != -32767)
    rf_mask = (rf_ds["lon"] >= rf_w) & (rf_ds["lon"] <= rf_e) & (rf_ds["lat"] >= rf_s) & (rf_ds["lat"] <= rf_n)
    rf_poi_mean = rf_clean.where(rf_mask).mean(dim=["y", "x"], skipna=True)
    rainforest_poi_wl = rf_poi_mean["wavelength_3d"].values / 1000.0
    rainforest_poi_refl = rf_poi_mean.values

    # Apply optional inverse correction to Galileo SSI
    if correction_factor is not None:
        Fig3["reflectance"] = Fig3["reflectance"] * correction_factor

    # PLOT
    fig, ax = plt.subplots(figsize=(7.08, 5), dpi=600)

    # PACE
    ax.plot(atacama_poi_wl, atacama_poi_refl, color="#e09f5f", linestyle="-", linewidth=1.5, alpha=0.7,)
    ax.plot(argentina_poi_wl, argentina_poi_refl, color="#a0a5ab", linestyle="-", linewidth=1.5, alpha=0.7,)
    ax.plot(rainforest_poi_wl, rainforest_poi_refl, color="#2e8b57", linestyle="-", linewidth=1.5, alpha=0.7,)

    # Galileo
    galileo_colors = {"a": "#e09f5f", "b": "#a0a5ab", "c": "#2e8b57"}
    galileo_lines = {"a": "--", "b": "--", "c": "--"}
    galileo_markers = {"a": "o", "b": "s", "c": "D"}

    for label, group in Fig3.groupby("dataset"):
        color = galileo_colors.get(label, "#888888")
        lines = galileo_lines.get(label, "--")
        markers = galileo_markers.get(label, "o")

        ax.plot(
            group["wl"], group["reflectance"],
            color=color, alpha=1, linestyle=lines, linewidth=1.75,
            marker=markers, markersize=4, markerfacecolor='white', markeredgewidth=1.5,
        )

    FONTSIZE = 12

    ax.tick_params(
        axis="both",
        which="major",
        labelsize=10,
    )

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='black', linestyle='--', label='Galileo SSI'),
        Line2D([0], [0], color='black', linestyle='-', label='PACE OCI')
    ]

    #ax.set_title("Earth Spectral Signatures: Modern PACE OCI vs. Galileo (1990)", fontsize=FONTSIZE, fontweight="bold", pad=12)
    ax.set_xlabel("Wavelength (μm)", fontsize=FONTSIZE)
    ax.set_ylabel("Reflectance", fontsize=FONTSIZE)
    ax.set_xlim(0.5, 1.050)
    ax.set_ylim(0, 0.45)
    ax.legend(frameon=True, fontsize=10, loc="lower right", ncol=1, handles=legend_elements)

    #plt.savefig("fig2.svg", format="svg", bbox_inches="tight")
    plt.tight_layout()
    plt.show()