"""Dictionary for handling PACE"""
#NOTE: spex_granules much narrower/don't overlap temporally with OCI & HARP2 (which are both a pretty close 1:1 match for FOV, and cover the same dates)

# PACE OCI SFREFL NRT (L2, V3.1) 5/10/26-6/20/26 & cloud_cover=(0,20) for granules (when possible); same for TOA (L1C, V3)
oci_granules = [
    {"name": "F3_B", "granule": "*20260628T172425*"}, # Region over Argentina, Galileo comparison
    {"name": "Pyramid", "granule": "*20251020T210415*"}, # Pyramid Lake cyanobacteria 
    {"name": "US_W", "granule": "*20260604T201346*"}, # GSL (UT), Craters (ID) 20251020T210415
    {"name": "US_W2", "granule": "*20260620T193954*"}, # Zuni Lave Fields (NM), Roger's Dry Lake (CA)
    {"name": "US_SE", "granule": "*20260514T175243*"}, # Everglades, swamps, temperate forests
    {"name": "NJ_blooms", "granule": "*20260604T165708*"}, # Coccolithophore/phytoplankton blooms (NJ, USA)
    {"name": "AU_W", "granule": "*20260618T050854*"}, # Shark Bay Stromatolites (AU)
    {"name": "AU_E", "granule": "*20260529T032329*"}, # Great Barrier Reef (NE Australia)
    {"name": "Mongolia", "granule": "*20260530T041317*"}, # Xilingol grasslands/steppes (Mongolia)
    {"name": "Amazon", "granule": "*20250827T170528*"}, # Amazon rainforest (Brazil; 0-20% cloud, no image available)
    {"name": "Atacama", "granule": "*20260511T173132*"}, # Atacame Desert (Chile)
    {"name": "Sahara", "granule": "*20260619T104836*"}, # Sahara Desert (Libya)
    {"name": "PA_Ocean1", "granule": "*20260613T203141*"}, # Pacific Ocean (near US W coast)
    {"name": "PA_Ocean2", "granule": "*20260613T233819*"}, # Pacific Ocean (open, southern hemisphere, 40% clouds)
    {"name": "Antarctica", "granule": "*20251031T152523*"}, # Antarctic granules ain't great in general
]

# PACE OCI (Sagan/Galileo comparison)
galileo_granules = [
    {"name": "PA_Ocean2", "granule": "*20260613T233819*"}, # Same open ocean as above, use for radiance
]
    # Unused granules for Galileo comparison...difficulties determining cloud coverage from Earthdata
    # {"name": "F1", "granule": "*20251208T020420*"}, # Fig.1 (Pacific Ocean cloud-free pixels) # PACE_OCI.20251208T020420.L1C.V3.5km.nc
    # {"name": "F3_AB", "granule": "*20251207T173742*"}, # Fig.3 Regions A & B (~cloud free)
    # {"name": "F3_C", "granule": "*20251207T174742*"}, # Fig.3 Region C (not cloud-free)

# PACE HARP2 V4, same granules as OCI (date/location matched)
harp_granules = [
    {"name": "US_W", "granule": "*20260604T201346*"}, # GSL, Salt Flats, Craters of the Moon, Roger's Dry Lake
    {"name": "US_W2", "granule": "*20260620T193954*"}, # Zuni Lava Fields
    {"name": "NJ_blooms", "granule": "*20260604T165708*"}, # Coccolithophores & Algal Bloom
    {"name": "US_SE", "granule": "*20260514T175243*"}, # Everglads, Okefenokee
    {"name": "Mongolia", "granule": "*20260530T041817*"}, # Xilingol Grasslands
    {"name": "Amazon", "granule": "*20250827T170528*"}, # Brazil
    {"name": "AU_E", "granule": "*20260529T032329*"}, # Great Barrier Reef
    {"name": "AU_W", "granule": "*20260618T050854*"}, # Stromatolites
    {"name": "Atacama", "granule": "*20260511T173132*"}, # Chile
    {"name": "Sahara", "granule": "*20260619T104836*"}, # Libya
    {"name": "PA_Ocean1", "granule": "*20260613T203141*"}, # near US W coast
    {"name": "PA_Ocean2", "granule": "*20260613T233819*"}, # open ocean, southern hemisphere (lots of surround cloud cover around poi)
    {"name": "Antarctica", "granule": "*20251031T152519*"}, # same polar issues
]

# PACE SPEXone V4, granules selected separately due to temporal differences from OCI/HARP2 (biome-specific, cloud cover unknown for many)
spex_granules = [
    {"name": "Craters", "granule": "*20240928T200026*"}, # Craters of the Moon (ID), also contains GSL
    # {"name": "Zuni", "granule": "*20260503T194104*"}, # Zuni Lave Fields (NM), # light clouds, ignore / find better granule
    {"name": "Rogers_Lake", "granule": "*20260420T202419*"}, # Roger's Dry Lake (CA) 
    {"name": "GSL", "granule": "*20260422T195547*"}, # Great Salt Lake (UT)
    {"name": "Salt_Flats", "granule": "*20260629T195727*"}, # Bonneville Salt Flats (UT), potentially some clouds
    {"name": "Coccolithophores", "granule": "*20260530T171950*"}, # (NJ) MAYBE contains coccolithophores cloud free 
    {"name": "Algal_Bloom", "granule": "*20260604T165708*"}, # (NJ) MAYBE contains algal bloom cloud free 
    {"name": "Everglades", "granule": "*20251227T180258*"}, # Everglades (FL)
    {"name": "Okefenokee", "granule": "*20260503T180245*"}, # Okefenokee Swamp (GA), minimal cloud cover
    {"name": "Amazon", "granule": "*20260502T171753*"}, # Amazon Rainforest (Brazil), cloudy / can't find alternative (ignore)
    {"name": "Mongolia", "granule": "*20260508T043845*"}, # Xilingol grasslands/steppes (Mongolia), #TODO: check for edge pixels
    {"name": "AU_W", "granule": "*20260120T053317*"}, # Shark Bay Stromatolites (NW AU)
    {"name": "AU_E", "granule": "*20251213T031039*"}, # Great Barrier Reef (NE AU), minimal clouds
    {"name": "Atacama", "granule": "*20260620T174635*"}, # Atacama Desert (Chile)
    {"name": "Sahara", "granule": "*20260511T111315*"}, # Sahara Desert (Libya)
    {"name": "PA_Ocean2", "granule": "*20260604T232024*"}, # Pacific Ocean (open, southern hemisphere), #TODO: change POI to avoid clouds
    {"name": "Antarctica", "granule": "*20260227T110237*"}, # Antarctic granules ain't great in general
]


# Sub-regions of interest from granules
biomes = {
    "US_W": [
        {
            "sub_name": "GSL",
            "class": "Microbial",
            "bbox": [-113, 40.6, -112, 41.8],  # [West, South, East, North]
            "poi": [-112.7, 41.1, -112.6, 41.35], # ~1:1 North:South arm of the GSL
            "color": "#ac0048",
        },
        {
            "sub_name": "GSL North Arm",
            "class": "Microbial",
            "bbox": [-113, 40.6, -112, 41.8],
            "poi": [-112.8, 41.24, -112.6, 41.4],
            "color": "#d67ea2",
        },
        {
            "sub_name": "GSL South Arm",
            "class": "Microbial",
            "bbox": [-113, 40.6, -112, 41.8],
            "poi": [-112.5, 40.95, -112.3, 41.15],
            "color": "#6dbfa4",
        },
        {
            "sub_name": "Wasatch Forest",
            "class": "Vegetation",
            "bbox": [-112, 40.0, -111, 42.0],
            "poi": [-111.85, 40.8, -111.6, 41.0], 
            "color": "#93E795",
        },
        {
            "sub_name": "Salt Flats", # Bonneville Salt Flats, W of GSL (UT)
            "class": "Abiotic",
            "bbox": [-114.2, 39.6, -113, 41.5],
            "poi": [-113.65, 39.9, -113.4, 40.1],
            "color": "#ACB173",
        },
        {
            "sub_name": "Craters of the Moon", # Craters of the Moon, ID, USA
            "class": "Abiotic",
            "bbox": [-113.8, 43, -113.3, 43.6],
            "poi": [-113.5, 43.225, -113.4, 43.35],
            "color": "#060000", #000000
        },
        {
            "sub_name": "Roger's Dry Lake", # Roger's Dry Lake, California (calcite-rich lake bed)
            "class": "Abiotic",
            "bbox": [-118.5, 34.5, -117.5, 35.5], # already very small
            "poi": [-117.9, 34.8, -117.75, 35],
            "color": "#DCE2C7",
        },
        {
            "sub_name": "Pyramid Lake (Summer)", # Pyramid Lake, NV, USA
            "class": "Microbial",
            "bbox": [-120, 39.5, -119, 40.5],
            "poi": [-119.65, 40, -119.5, 40.1],
            "color": "#89AD9C",
        }
    ],
    "Pyramid": [
        {
            "sub_name": "Pyramid Lake (Fall)", # Pyramid Lake, NV, USA (Peak cyanobacterial bloom season)
            "class": "Microbial",
            "bbox": [-120, 39.5, -119, 40.5],
            "poi": [-119.65, 40, -119.5, 40.1],
            "color": "#206345",
        }
    ],
    "US_W2": [
        {
            "sub_name": "Zuni Lava Fields", # Zuni-Bandera volcanic field, NM, USA
            "class": "Abiotic",
            "bbox": [-108.15, 34.7, -107.8, 35.1],
            "poi": [-107.95, 34.86, -107.9, 34.94],
            "color": "#3F3F3F",
        }
    ],
    "NJ_blooms": [
        {
            "sub_name": "Coccolithophores",
            "class": "Microbial",
            "bbox": [-74.5, 38, -67, 42], # contains some land; -74.5, 38, -71, 41
            "poi": [-73.4, 39.2, -73.2, 39.6],
            "color": "#ace1cf",
        },
        {
            "sub_name": "Algal Bloom",
            "class": "Microbial",
            "bbox": [-74.5, 38, -67, 42], # contains some land; -71, 40, -67, 42.5
            "poi": [-68.3, 41.6, -68.1, 41.7],
            "color": "#05783F",
        },
    ],
    "US_SE": [
        {
            "sub_name": "Everglades",
            "class": "Vegetation",
            "bbox": [-81, 25.1, -80.5, 25.7],
            "poi": [-80.9, 25.4, -80.7, 25.6],
            "color": "#2C7A1A",
        },
        {
            "sub_name": "Okefenokee",
            "class": "Vegetation",
            "bbox": [-82.7, 30.5, -82.1, 31],
            "poi": [-82.6, 30.8, -82.3, 31],
            "color": "#59E939",
        },
    ],
    "Mongolia": [
        {
            "sub_name": "Grasslands", # Xilin Gol Grasslands
            "class": "Vegetation",
            "bbox": [115.1, 43.0, 116.5, 44.2],
            "poi": [115.8, 43.4, 116.2, 43.8],
            "color": "#AFC16A",
        }
    ],
    "Amazon": [
        {
            "sub_name": "Rainforest", # Amazon rainforest near Serra da Cutia National Park (Brazil)
            "class": "Vegetation",
            "bbox": [-65.0, -10.0, -60.0, -5.0],
            "poi": [-62.5, -8.4, -62.0, -8.2],
            "color": "#042D00",
        }
    ],
    "AU_E": [
        {
            "sub_name": "Coral Reef", # Great Barrier Reef, NE Australia
            "class": "Microbial",
            "bbox": [149, -20, 150, -19.5], # 145, -20, 150, -16
            "poi": [149.35, -19.8, 149.5, -19.7], # 148, -19.5, 149, -19
            "color": "#6D7BFF",
        },
    ],
    "AU_W": [
        {
            "sub_name": "Stromatolites", # Shark Bay Stromatolites, Western Australia #TODO: compare with GSL stromatolites
            "class": "Microbial",
            "bbox": [113.2, -26.7, 114.3, -25.3], 
            "poi": [114, -26.4, 114.2, -26.2],
            "color": "#B6BDFF",
        },
    ],
    "Atacama": [
        {
            "sub_name": "Atacama", # Atacama Desert, Chile
            "class": "Abiotic",
            "bbox": [-70.0, -24.5, -68.5, -22.0],
            "poi": [-69.4, -23.7, -69.2, -23.5], # -69.5, -24, -69, -23.5 (larger area)
            "color": "#677341",
        }
    ],
    "Sahara": [
        {
            "sub_name": "Sahara", # Libyan desert
            "class": "Abiotic",
            "bbox": [21.0, 24.0, 24.0, 27.0],
            "poi": [21.4, 24.7, 21.6, 25], # 21.2, 24.5, 21.8, 25.5 (larger area)
            "color": "#900000",
        }
    ],
    "PA_Ocean1": [
        {
            "sub_name": "CA Upwelling",
            "class": "Abiotic",
            "bbox": [-136, 35, -125, 45],
            "poi": [-127, 39.8, -125.8, 40.8], # focuses on bloom, also less active areas in granule
            "color": "#247A97",
        },
        {
            "sub_name": "CA Opean Ocean",
            "class": "Abiotic",
            "bbox": [-136, 35, -125, 45],
            "poi": [-133.5, 37.5, -132, 38.5], # focuses on open ocean
            "color": "#288CAD",
        },
    ],
    "PA_Ocean2": [
        {
            "sub_name": "Open Ocean",
            "class": "Abiotic",
            "bbox": [-140.0, -35.0, -135.0, -30.0],
            "poi": [-155.5, -6, -154, -5],
            "color": "#1F5E75",
        }
    ],
    "Antarctica": [
        {
            "sub_name": "Antarctica",
            "class": "Abiotic",
            "bbox": [60, -87, 70, -80],
            "poi": [64.5, -86, 65, -85], # 64, -86, 65, -85 (larger area)
            "color": "#B9E3F1",
        }
    ],
}


# GSL Regions of interest for case study analysis
GSL_biomes = [
    {
        "name": "(A) GSL North Arm",
        "bbox": [-113, 40.6, -112, 41.8],
        "poi": [-112.8, 41.24, -112.6, 41.4],
        "color": "#d67ea2",
    },
    {
        "name": "(B) GSL South Arm",
        "bbox": [-113, 40.6, -112, 41.8],
        "poi": [-112.5, 40.95, -112.3, 41.15],
        "color": "#6dbfa4",
    },
    {
        "name": "(C) Wasatch Forest",
        "bbox": [-112, 40.0, -111, 42.0],
        "poi": [-111.85, 40.8, -111.6, 41.0], 
        "color": "#015C02",
    },
    {
        "name": "(D) Salt Flats",
        "bbox": [-114.2, 39.6, -113, 41.5],
        "poi": [-113.65, 39.9, -113.4, 40.1], 
        "color": "#CDCDC5",
    },
    {
        "name": "Bounding Box",
        "bbox": [-113.25, 40.2, -111.5, 42],
        "poi": [-113.25, 40.2, -111.5, 42],
        "color": "#000000",
    },
]


# Biosignatures of interest:

# Gas dictionary for calculating volume mixing ratios (VMRs; in um for Sagan comparison)
gases_config = {
    "O2": {
        "center": 0.762,     # O2-A band absorption core
        "L": 0.745,     # Near-continuum left window
        "R": 0.780,    # Near-continuum right window
        "cross_section": 5.9e-26, # cm^2/molecule convolved to 5nm
    },
    "H2O": {
        "center": 0.940,     # Strong Water vapor absorption
        "L": 0.880,     # Clean window left
        "R": 1.010,    # Clean window right
        "cross_section": 1.6e-24, # cm^2/molecule
    },
    "O3": {
        "center": 0.600,     # Broad Chappuis band peak
        "L": 0.510,     # Continuum baseline left
        "R": 0.680,    # Continuum baseline right
        "cross_section": 5.1e-21, # cm^2/molecule
    },
}

# Absorption band dictionary for pigment biosignatures of interest (GSL case study; nm)
bands_config = {
    # Gas Biosignatures / Atmospheric Absorption
    # "O2-B": {
    #     "center": 688,
    #     "band": (685, 692),
    #     "L": (675, 682),
    #     "R": (695, 702),
    # },
    "O2-A": {
        "center": 760,
        "band": (758, 770),
        "L": (748, 752),
        "R": (772, 776),
    },
    # "H2Oa": {
    #     "center": 720,
    #     "band": (715, 730),
    #     "L": (700, 710),
    #     "R": (735, 745),
    # },
    "H2Ob": {
        "center": 820,
        "band": (812, 830),
        "L": (800, 806),
        "R": (838, 846),
    },
    # "H2Oc": {
    #     "center": 940,
    #     "band": (925, 960),
    #     "L": (900, 915),
    #     "R": (970, 985),
    # },
    # Microbial Pigments (in vivo shifts from reference features based on observed spectra)
    "Bacterioruberin": {
        "center": 540,
        "band": (520, 560),
        "L": (490, 510),
        "R": (570, 590),
    },
    "Phycocyanin": {
        "center": 620,
        "band": (600, 640),
        "L": (570, 590),
        "R": (650, 670),
    },
    "β-carotene": {
        "center": 500,
        "band": (485, 515),
        "L": (440, 460),
        "R": (530, 550),
    },
    "Chlorophyll A": {
        "center": 680,
        "band": (665, 690),
        "L": (610, 630),
        "R": (710, 730),
    },
}