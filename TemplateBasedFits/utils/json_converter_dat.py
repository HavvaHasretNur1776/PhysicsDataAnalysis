import json
import numpy as np
import os

# ---- Input JSON files ----
json_files = [
    "histos_2016_up.json",
    "histos_2016_down.json",
    "histos_2017_up.json",
    "histos_2017_down.json",
    "histos_2018_up.json",
    "histos_2018_down.json"
]

# ---- Map JSON keys to your fitter names ----
'''component_map = {
    "b0_dstmunu": "DstMu",
    "b0_dstdx": "DstDX",
    "b0_dststmu": "DstDst",
}'''

# ---- Load all JSON files ----
all_data = []
for file in json_files:
    with open(file) as f:
        all_data.append(json.load(f))

# Get components automatically from first year
components = all_data[0].keys()

# ---- Output directory ----
output_dir = "./"
os.makedirs(output_dir, exist_ok=True)

# ---- Loop over components ----
#for json_key, comp_name in component_map.items():
for comp in components:

    # assume all years have same wbin structure
#    wbins = all_data[0][json_key].keys()
    wbins = all_data[0][comp].keys()

    for wbin in wbins:

        # Sum histograms across years
        combined = None

        for year_data in all_data:
#            values = np.array(year_data[json_key][wbin]["Value"], dtype=float)
            values = np.array(year_data[comp][wbin]["Value"], dtype=float)

            if combined is None:
                combined = values
            else:
                combined += values

        # Write output
#        filename = f"{comp_name}_wbin{wbin}.dat"
        filename = f"{comp}_wbin{wbin}.dat"
        filepath = os.path.join(output_dir, filename)

        with open(filepath, "w") as f_out:
            for v in combined:
                f_out.write(f"{v}\n")

        print(f"Written {filepath}")