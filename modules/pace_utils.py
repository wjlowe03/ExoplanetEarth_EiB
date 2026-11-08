"""Process/Plot PACE Instrumental Data"""

import os
import numpy as np
import xarray as xr
import pandas as pd
import earthaccess # type: ignore
import pickle
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from matplotlib.patches import Patch



# PROCESS: OCI L2 SFREFL (cloud/ice masking available)
def _process_sfrefl(path):
    """Processes OCI Surface Reflectance into a standardized dataset layout."""
    dt = xr.open_datatree(path)
    ds = xr.merge(
        (
            dt.ds,
            dt["geophysical_data"].ds[["rhos", "l2_flags"]],
            dt["sensor_band_parameters"].coords,
            dt["navigation_data"].ds.set_coords(
                ("longitude", "latitude")
            ).coords,
        )
    )

    # # CHECK
    # print(ds["l2_flags"])
    # print(ds["l2_flags"].attrs)

    # Cloud masking to protect terrestrial surface averages
    cloud = get_l2_flag(ds["l2_flags"], "CLDICE")
    rhos_masked = ds["rhos"].where(~cloud)

    # cld_mask = ~(ds["l2_flags"].cf == "CLDICE")
    # rhos_masked = ds["rhos"].where(cld_mask)

    target = rhos_masked.to_dataset(name="data_var")
    return target.rename(
        {
            "latitude": "lat",
            "longitude": "lon",
            "number_of_lines": "y",
            "pixels_per_line": "x",
        }
    )


# PROCESS: OCI L1C TOA Reflectance
def _process_reflectance(path, nadir_idx=0): 
    """Processes L1C file into a TOA reflectance dataset.""" 
    # NOTE: If files come back all NA, change to nadir_idx=1

    dt = xr.open_datatree(path)

    with xr.open_dataset(path) as root:
        d_earth_sun = float(root.attrs.get("sun_earth_distance", 1.0))

    obs_view = dt["observation_data"].ds.isel(number_of_views=nadir_idx)
    geo_view = dt["geolocation_data"].ds.isel(number_of_views=nadir_idx)
    view_angle_profile = dt["sensor_views_bands"].ds.isel(number_of_views=nadir_idx)

    toa_radiance = obs_view["i"]  # (bins_along_track, bins_across_track, intensity_bands_per_view)
    toa_radiance = toa_radiance.where(toa_radiance >= 0, np.nan)
    solar_zenith_deg = geo_view["solar_zenith_angle"]  # (bins_along_track, bins_across_track)
    f0 = view_angle_profile["intensity_f0"]  # (intensity_bands_per_view)
    dim_axis_name = "intensity_bands_per_view"
    wv_var_name = "intensity_wavelength"

    wavelengths = view_angle_profile[wv_var_name]

    # Compute TOA Reflectance: (pi * L * d^2) / (F0 * cos(SZA))
    solar_zenith_rad = np.radians(solar_zenith_deg)
    cos_sza = np.cos(solar_zenith_rad)
    cos_sza = cos_sza.where(cos_sza > 0)  # drop zero values to avoid errors
    
    numerator = np.pi * toa_radiance * (d_earth_sun**2)
    denominator = f0 * cos_sza
    toa_reflectance = numerator / denominator

    toa_ds = toa_reflectance.to_dataset(name="rhos")

    # Remerge/rename everything
    ds = xr.merge((toa_ds, geo_view[["latitude", "longitude"]]))

    ds = ds.rename(
        {
            "latitude": "lat",
            "longitude": "lon",
            "bins_along_track": "y",
            "bins_across_track": "x",
            dim_axis_name: "wavelength_3d",
        }
    )

    ds = ds.assign_coords(wavelength_3d=wavelengths.values)
    ds = ds.set_coords(["lat", "lon"])

    target = ds["rhos"].to_dataset(name="data_var")  # isolate rhos but keep lat/lon/wavelength_3d coords intact

    return target


# PROCESS: OCI L1C TOA Radiance
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


# PROCESS: SPEXone L1C SCI Polarization
def _process_spex(path):
    view = xr.open_dataset(path, group="sensor_views_bands").squeeze()
    geo = xr.open_dataset(path, group="geolocation_data")
    obs = xr.open_dataset(path, group="observation_data").squeeze()

    view_angles = view["sensor_view_angle"].values
    pol_wavelengths = view["polarization_wavelength"].values
    int_wavelengths = view["intensity_wavelength"].values

    # Explicitly map the coordinate to the correct underlying dimension (account for different dimensions)
    ds_intensity = xr.Dataset(
        {"i": obs["i"]},
        coords={
            "number_of_views": ("number_of_views", view_angles),
            "intensity_wavelength": (
                ("number_of_views", "intensity_bands_per_view"), 
                int_wavelengths,
            ),
        },
    )

    ds_pol = xr.Dataset(
        {"q": obs["q"], "u": obs["u"], "dolp": obs["dolp"]},
        coords={
            "number_of_views": ("number_of_views", view_angles),
            "polarization_wavelength": (
                ("number_of_views", "polarization_bands_per_view"),
                pol_wavelengths,
            ),
        },
    )

    i_interpolated_list = []
    for v in range(len(view_angles)):
        v_int = ds_intensity.isel(number_of_views=v)
        v_pol = ds_pol.isel(number_of_views=v)

        # Swap index dimension for wavelength dimension (xarray interpolation)
        v_int_swapped = v_int["i"].swap_dims({"intensity_bands_per_view": "intensity_wavelength"})
        v_i_interp = v_int_swapped.interp(intensity_wavelength=v_pol["polarization_wavelength"].values, method="linear",)
        v_i_interp = v_i_interp.rename({"intensity_wavelength": "polarization_bands_per_view"})
        i_interpolated_list.append(v_i_interp)

    # Recombine the interpolated intensity across all views, merge
    i_interp_all_views = xr.concat(i_interpolated_list, dim="number_of_views")
    ds = xr.merge([ds_pol, i_interp_all_views.to_dataset(name="i")])

    # Add Geolocation & Standardize Names
    ds = ds.assign_coords(
        view_angle=("number_of_views", view_angles),
        lat=geo["latitude"],
        lon=geo["longitude"],
    )

    return ds.rename(
        {
            "bins_along_track": "y",
            "bins_across_track": "x",
            "number_of_views": "view",
        }
    )


# PROCESS: HARP2 L1C SCI Polarization
def _process_harp(path):

    view = xr.open_dataset(path, group="sensor_views_bands").squeeze()
    geo  = xr.open_dataset(path, group="geolocation_data")
    obs  = xr.open_dataset(path, group="observation_data").squeeze()

    ds = obs.rename(
        {
            "bins_along_track": "y",
            "bins_across_track": "x",
            "number_of_views": "view",
        }
    )

    ds = ds.assign_coords(
        view=np.arange(ds.sizes["view"]),
        view_angle=("view", view["sensor_view_angle"].values),
        wavelength=("view", view["intensity_wavelength"].values),
        lat=(("y","x"), geo["latitude"].values),
        lon=(("y","x"), geo["longitude"].values),
    )

    return ds


# Instrument registry for extracting spectra from earthaccess, cache to .pkl
_INSTRUMENT_REGISTRY = {

    "PACE_OCI_L2_SFREFL": {
        "processor": _process_sfrefl,
        "granule_key": "granule",
        "mean_variable": "data_var",
    },

    "PACE_OCI_L2_SFREFL_NRT": {
        "processor": _process_sfrefl,
        "granule_key": "granule",
        "mean_variable": "data_var",
    },

    "PACE_OCI_L1C_SCI": {
        "processor": _process_reflectance,
        "granule_key": "granule",
        "mean_variable": "data_var",
    },

    # "PACE_OCI_L1C_SCI": {
    #     "processor": _process_radiance,
    #     "granule_key": "granule",
    #     "mean_variable": "data_var",
    # },

    "PACE_SPEXONE_L1C_SCI": {
        "processor": _process_spex,
        "granule_key": "granule",
        "mean_variable": "dolp",
    },

    "PACE_HARP2_L1C_SCI": {
        "processor": _process_harp,
        "granule_key": "granule",
        "mean_variable": "dolp",
    },
}


# STREAM: Extract PACE spectra from earthaccess based on instrument/granule/processing function
def extract_spectra(short_name, granule_dict, cache_dir=None, **kwargs):
    
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
        # Process the remote granule
        ds = instrument["processor"](paths[0], **kwargs)
        mean_var = instrument["mean_variable"]

        ds.load()
        ds.close()

        # Calculate summary statistics from the fully loaded dataset
        spectra_results = {
            "granule_mean": ds[mean_var].mean(
                dim=["y", "x"], skipna=True
            ),
            "subsample_mean": None,
            "biome_mean": None,
            "ds_reference": ds,
        }

        print(f"Processing completed for '{item_name}'.")

        # Cache the fully materialized results
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


# Call on cached files
def get_cached_files(cache_dir, biomes):

    valid_files = []
    total_profiles = 0

    if not os.path.exists(cache_dir):
        return valid_files, total_profiles

    for filename in sorted(os.listdir(cache_dir)):

        if not filename.endswith("_cache.pkl"):
            continue

        base = filename.replace("_cache.pkl", "")

        if "_PACE_" in base:
            parent_name = base.split("_PACE_")[0]
        else:
            parent_name = base

        regions = biomes.get(parent_name, [])
        total_profiles += len(regions)
        valid_files.append((filename, parent_name))

    return valid_files, total_profiles


def get_l2_flag(flags, flag_name):

    names = flags.attrs["flag_meanings"].split()
    masks = flags.attrs["flag_masks"]
    lookup = dict(zip(names, masks))

    if flag_name not in lookup:
        raise ValueError(f"Unknown flag '{flag_name}'")

    return (flags & lookup[flag_name]) != 0


def pkl_to_csv(data_dir, file_name, output_csv):
    """Convert a PACE .pkl cache file to .csv."""

    cache_file = os.path.join(data_dir, file_name)

    print("Loading pickle file...")
    with open(cache_file, "rb") as f:
        cache_data = pickle.load(f)

    ds = cache_data["ds_reference"]

    # Flatten & remove empty rows
    print("Reshaping data...")
    clean_data = ds["data_var"].where(ds["data_var"] != -32767)
    ds_flat = clean_data.stack(pixel=("y", "x"))
    ds_flat = ds_flat.dropna(dim="pixel", how="all")

    # Extract matrix & transpose values
    lats = ds_flat["lat"].values
    lons = ds_flat["lon"].values
    wavelengths = ds_flat["wavelength_3d"].values
    reflectance_matrix = ds_flat.values.T

    # Build dataframe
    print("Building dataframe...")
    df_final = pd.DataFrame(reflectance_matrix, columns=wavelengths)
    df_final.insert(0, "lon", lons)
    df_final.insert(0, "lat", lats)

    df_final.to_csv(output_csv, index=False)

    print(f"Done! Saved to: {output_csv}")
    print(f"Shape: {df_final.shape}")

    return df_final


# # Plot OCI spectra/maps
def plot_biome_spectra(
    valid_files, cache_dir, biomes, target_sub_biomes=None, spatial_level="bbox", xlim=None,
):

    if isinstance(spatial_level, str):
        spatial_levels = [spatial_level]
    else:
        spatial_levels = spatial_level

    # Set allowed spatial levels
    allowed_levels = {"granule", "bbox", "poi"}
    if not set(spatial_levels).issubset(allowed_levels):
        raise ValueError(f"spatial_level must be one or more of {allowed_levels}")
    
    # Set line styles for plotting
    line_styles = {"granule": ":", "bbox": "--", "poi": "-"}

    # Plot
    fig, ax = plt.subplots(figsize=(6, 4), dpi=600) #13, 8

    plotted_regions = set()
    n_plotted = 0

    biome_spectra_dict = {}

    # Loop through the cached granules
    for filename, parent_name in valid_files:
        cache_file_path = os.path.join(cache_dir, filename)

        with open(cache_file_path, "rb") as f:
            data = pickle.load(f)

        ds = data["ds_reference"]

        for region in biomes.get(parent_name, []):
            sub_name = region["sub_name"]

            if (
                target_sub_biomes is not None
                and sub_name not in target_sub_biomes
            ):
                continue

            for level in spatial_levels:

                # Determine "region" keys (for granule vs. poi)
                if level == "granule":
                    region_key = (filename, level)
                else:
                    if level not in region:
                        continue
                    if not region[level]:
                        continue
                    region_key = (filename, level, tuple(region[level]))

                # Avoid duplicate granule/bbox plots
                if region_key in plotted_regions:
                    continue

                plotted_regions.add(region_key)

                # Calculate mean spectrum for various spatial levels
                if level == "granule":
                    valid_mask = ~np.isnan(ds["data_var"].isel(wavelength_3d=0))
                    pixel_count = int(valid_mask.sum())
                    mean_profile = (ds["data_var"].mean(dim=["y", "x"], skipna=True))
                    label = (f"{parent_name} (granule)")

                elif level == "bbox":
                    w, s, e, n = region["bbox"]
                    spatial_mask = ((ds["lon"] >= w) & (ds["lon"] <= e) & (ds["lat"] >= s) & (ds["lat"] <= n))
                    valid_mask = (spatial_mask & ~np.isnan(ds["data_var"].isel(wavelength_3d=0)))
                    pixel_count = int(valid_mask.sum())
                    mean_profile = (ds["data_var"].where(spatial_mask).mean(dim=["y", "x"], skipna=True))
                    label = (f"{sub_name} (bbox)")

                elif level == "poi":
                    w, s, e, n = region["poi"]
                    spatial_mask = ((ds["lon"] >= w) & (ds["lon"] <= e) & (ds["lat"] >= s) & (ds["lat"] <= n))
                    valid_mask = (spatial_mask & ~np.isnan(ds["data_var"].isel(wavelength_3d=0)))
                    pixel_count = int(valid_mask.sum())
                    mean_profile = (ds["data_var"].where(spatial_mask).mean(dim=["y", "x"], skipna=True))
                    label = (f"{sub_name} (poi)")

                # Validate spectrum
                valid_pixels = np.count_nonzero(~np.isnan(mean_profile.values))

                if pixel_count == 0:
                    print(f"{label}: no contributing spatial pixels")
                    continue

                color = region.get("color", None)
                print("Plotting:", label, "| spatial pixels:", pixel_count)

                ax.plot(
                    mean_profile["wavelength_3d"].values,
                    mean_profile.values,
                    color=color,
                    linestyle=line_styles[level],
                    linewidth=1.5,
                    label=label,
                )

                n_plotted += 1

    ax.tick_params(axis='both', labelsize=8)
    ax.set_xlabel("Wavelength (nm)", fontsize=12)
    ax.set_ylabel("Reflectance", fontsize=12)

    if xlim is not None:
        ax.set_xlim(*xlim)

    if n_plotted > 0:
        ax.legend(title="Biome", ncol=2, loc="upper left", fontsize=6)

    else:
        print("No spectra plotted.")

    plt.tight_layout()
    plt.show()

    return biome_spectra_dict


# PLOT: OCI SF and TOA reflectances by spatial reference
def plot_biome_maps(
    valid_files, cache_dir, biomes, target_wl, target_sub_biome, bbox_override=None, cmap="viridis",
):

    for filename, parent_name in valid_files:
        cache_file_path = os.path.join(cache_dir, filename)

        with open(cache_file_path, "rb") as f:
            data = pickle.load(f)

        ds = data["ds_reference"]

        # Find requested sub-biome
        for region in biomes.get(parent_name, []):

            if region["sub_name"] != target_sub_biome:
                continue

            # Select nearest wavelength
            image = ds["data_var"].sel(wavelength_3d=target_wl, method="nearest")
            actual_wavelength = float(image["wavelength_3d"].values)

            # Plot
            fig, ax = plt.subplots(figsize=(7, 6))

            # Crop to region of interest
            if bbox_override is not None:
                w, s, e, n = bbox_override
                spatial_mask = ((ds["lon"] >= w) & (ds["lon"] <= e) & (ds["lat"] >= s)& (ds["lat"] <= n))

                if not np.any(spatial_mask.values & ~np.isnan(image.values)):
                    print(f"No valid pixels found within custom bbox.")
                    plt.close(fig)
                    return

                plot_image = image.where(spatial_mask)

                ax.set_xlim(w, e)
                ax.set_ylim(s, n)
                ax.set_aspect("equal", adjustable="box")
                title_region = "Custom ROI"

            # Plot full granule
            else:
                plot_image = image
                title_region = "Full Granule"

            # Plot
            pcm = ax.pcolormesh(
                ds["lon"],
                ds["lat"],
                plot_image,
                shading="auto",
                cmap=cmap,
            )

            plt.colorbar(pcm, ax=ax, label="Reflectance")

            ax.set_title(
                f"{parent_name} → {target_sub_biome}\n"
                f"{title_region}\n"
                f"{actual_wavelength:.1f} nm",
                fontsize=11,
                weight="bold",
            )

            ax.set_xlabel("Longitude")
            ax.set_ylabel("Latitude")
            ax.grid(True, linestyle=":", alpha=0.3)

            plt.tight_layout()
            plt.show()

            return

    print(f"Sub-biome '{target_sub_biome}' not found.")


def plot_biomes_csv(csv, biomes, filter_mode="poi"): 
    """Plot mean spectra for each biome and the whole granule.""" 
 
    df_master = csv
 
    dfs = {} 
    colors = {} 
    biome_names = [] 
 
    # Loop through biomes to catalog/splice/etc. 
    for b in biomes: 
        name = b["name"] 
        color = b["color"] 
        w, s, e, n = b[filter_mode] 
        spatial_mask = ( 
            (df_master["lon"] >= w) & (df_master["lon"] <= e) & 
            (df_master["lat"] >= s) & (df_master["lat"] <= n) 
        ) 
        df_slice = df_master[spatial_mask] 
        if len(df_slice) > 0: 
            biome_names.append(name) 
            dfs[name] = df_slice 
            colors[name] = color 
        else: 
            print(f"Warning: No overlapping data coordinates found for {name}") 
 
    # Append entire granule as a "biome" 
    biome_names.append("Whole Granule") 
    dfs["Whole Granule"] = df_master 
    colors["Whole Granule"] = "black" 

    # Plot
    fig, ax = plt.subplots(figsize=(3, 2), dpi=600)

    location_handles = []

    for biome, df in dfs.items():
        wavelengths = np.array(df.columns[2:], dtype=float)
        spectra = df.iloc[:, 2:].to_numpy(dtype=float)

        # Mean spectra across all pixels
        mean_spectrum = np.nanmean(spectra, axis=0)
        valid = ~np.isnan(mean_spectrum)

        if biome == "Whole Granule":
            line, = ax.plot(
                wavelengths[valid], mean_spectrum[valid],
                color=colors[biome], linestyle=":", linewidth=1.0
            )
        elif biome == "Bounding Box":
            line, = ax.plot(
                wavelengths[valid], mean_spectrum[valid],
                color=colors[biome], linestyle="--", linewidth=1.0
            )
        else:
            line, = ax.plot(
                wavelengths[valid], mean_spectrum[valid],
                color=colors[biome], linewidth=1.0
            )

            location_handles.append(
                Line2D([0], [0], color=colors[biome],
                    linewidth=1.0, label=biome)
            )

    # Custom handles for line styles

    linetype_handles = [
        Line2D([0], [0], color="black", linestyle=":",
            linewidth=1.0, label="Whole Granule"),
        Line2D([0], [0], color="black", linestyle="--",
            linewidth=1.0, label="Bounding Box")
    ]

    # Axis formatting
    ax.set_xlim(345, 900)
    ax.set_xlabel("Wavelength (nm)", fontsize=8)
    ax.set_ylabel("Reflectance", fontsize=8)
    ax.tick_params(axis="both", labelsize=8)

    # First legend: biome/location identifiers
    leg1 = ax.legend(
        handles=location_handles,
        loc="upper center",
        bbox_to_anchor=(0.25, -0.25),
        fontsize=6,
        frameon=False,
        ncol=1
    )
    ax.add_artist(leg1)

    # Second legend: line styles
    ax.legend(
        handles=linetype_handles,
        loc="upper center",
        bbox_to_anchor=(0.75, -0.25),
        fontsize=6,
        frameon=False,
        ncol=1
    )

    # Reserve space below the plot for legends
    fig.subplots_adjust(
        left=0.15,
        right=0.95,
        bottom=0.35,
        top=0.95
    )

    # fig.savefig("biomes.png", bbox_inches="tight", dpi=600)
    plt.show()