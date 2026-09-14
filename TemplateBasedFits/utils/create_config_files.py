import json
import os
import copy

# Load the template config (your bin1 version)
with open("config_bin1.json") as f:
    config_template = json.load(f)

# Example custom mmiss2 binning per bin (update as needed)
mmiss2_binning_per_bin = {
   1: [-4.0, 5.5],
    2: [-5.0, 6.0],
    3: [-3.0, 5.5],
    4: [-3.0, 5.0],
    5: [-3.0, 5.0],
    6: [-3.0, 5.0],
    7: [-2.8, 5.0],
    8: [-2.0, 5.0],
    9: [-2.0, 4.1],
    10: [-2.0, 3.0],
}

# Loop over bins 1 to 10 and create config files
for i in range(1, 11):
    config = copy.deepcopy(config_template)

    for comp in config.values():
        # Replace _bin1 with _bin{i} in file paths
        comp["fileNames"] = [fn.replace("bin1", f"bin{i}") for fn in comp["fileNames"]]
        # Replace histoname
        comp["histoname"] = f"mmiss2rec_brf_histo_20bin_wbin{i}"
        # Set binning
        comp["binning"] = [mmiss2_binning_per_bin[i]]

    # Save to new config file
    out_name = f"config_bin{i}.json"
    with open(out_name, "w") as out_f:
        json.dump(config, out_f, indent=2)

    print(f"Saved: {out_name}")
