import numpy as np
import matplotlib.pyplot as plt
import os 
import mplhep as hep
import json

# Set LHCb plotting style
hep.style.use("LHCb2")
# Define colors
colors = plt.cm.tab10.colors

# Output directory
#output_dir = "all_templates_plots"
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

#base_path = "merged_years"
handles = []
yields = {}
#components = ["nDstMu", "nDstTau", "nDstDX", "nb_dstst", "nbs_dstst", "nDststHigher","nWSPi","nWSMu"]#"WSPi","WSMu","WSMuWSPi","WSPi-WSMuWSPi"]
component_labels = {
    "b0_dstmunu": r"$B^0 \to D^{*+} \mu^- \nu_\mu$",
    "b0_dsttaunu": r"$B^0 \to D^{*+} \tau^- \nu_\tau$",
    "b_dstdx": r"$B^+ \to D^{*+}D(X_c \to \mu \nu X)$",
    "b_dststmu_heavy": r"$B \to D^{**}(Heavy) \mu^- \nu_\mu$",
    "WSMu": r"Combinatorial",
    "WSPi": r"$D^{*}$ Combinatorial",
    "b_dststmu": r"$B \to D^{**}(Light) \mu \nu$",
    "bs_dststmu": r"$B_s \to D_s^{**} \mu \nu$"
}

'''filenames = {
    'nDstMu' : 'b0_dstmunu_wbin' ,
    'nDstTau' : 'b0_dsttaunu_wbin' ,

    'nDstDX' : 'nDstDX_bin',
    'nb_dstst' : 'nb_dstst_bin',
    'nbs_dstst' : 'nbs_dstst_bin',
    'nDststHigher' : 'nDststHigher_bin',

    'nWSPi' :  'WSPi_bin',
    'nWSMu' :  'WSMu_bin',
    'nWSPi-WSMuWSPi' :  'WSPi-WSMuWSPi_bin'
}
'''
# Settings
nbins = 20
custom_mmiss2_bins = [
        [-4., 5.5],   # Bin 1
        [-5., 6.],    # Bin 2
        [-3., 5.5],    # Bin 3
        [-3., 5.],     # Bin 4
        [-3., 5.],     # Bin 5
        [-3., 5.],     # Bin 6
        [-2.8, 5.],   # Bin 7
        [-2., 5.], # Bin 8
        [-2., 4.1], # Bin 9
        [-2., 3.]  # Bin 10
    ]

'''custom_mmiss2_bins = [
        [-4., 10.5],   # Bin 1
        [-5., 8.],    # Bin 2
        [-3., 8.5],    # Bin 3
        [-3., 7.],     # Bin 4
        [-3., 7.],     # Bin 5
        [-3., 7.],     # Bin 6
        [-2.8, 7.],   # Bin 7
        [-2., 7.], # Bin 8
        [-2., 6.1], # Bin 9
        [-2., 10.]  # Bin 10
    ]'''

yield_name_map = {
    'b0_dstmunu': 'nDstMu',
    'b0_dsttaunu': 'nDstTau',
    'b_dstdx': 'nDstDX',
    'b_dststmu': 'nb_dstst',
    'b_dststmu_heavy': 'nDststHigher',
    'bs_dststmu': 'nbs_dstst',
    'WSPi': 'nWSPi',
    'WSMu': 'nWSMu'
}
components = ['b0_dstmunu', 'b0_dsttaunu', 
#                'b0_dstdx', 'b+_dstdx','b0_dststmu_heavy', 'b+_dststmu_heavy',
                'b_dstdx', 'b_dststmu_heavy','b_dststmu','bs_dststmu',
                #'b0_dststmu_10413', 'b0_dststmu_20413', 'b0_dststmu_415',
                #'b+_dststmu_10423', 'b+_dststmu_20423', 'b+_dststmu_425',
                #'bs_dststmu_10433', 'bs_dststmu_435',
                'WSPi','WSMu']

# Loop over bins 1 to 10
for bin_index in range(0, 10):
    print(f"BIN {bin_index}")
    plt.figure(figsize=(8, 6))
    handles = []
    for i, comp in enumerate(components):
        if comp in ['b0_dstmunu', 'b0_dsttaunu', 'WSPi','WSMu']:
            file_path = f"./{comp}_wbin{bin_index}.dat"
        else:
            file_path = f"grouped_templates/{comp}_wbin{bin_index}.dat"

#        bin_edges = np.linspace(custom_mmiss2_bins[bin_index-1][0], custom_mmiss2_bins[bin_index-1][1], 21)  
        bin_edges = np.linspace(custom_mmiss2_bins[bin_index][0], custom_mmiss2_bins[bin_index][1], 21)  

        hist_values = np.loadtxt(file_path)#
        
#        yields[comp] = float(np.sum(hist_values))
        yield_name = yield_name_map[comp]
        yields[yield_name] = float(np.sum(hist_values))
        # Normalize the histogram values
        hist_values = hist_values / np.sum(hist_values) if np.sum(hist_values) > 0 else hist_values
        #line, = plt.step(bin_edges[:-1], hist_values, where='post', label=label, color=color)
        line, = plt.step(bin_edges[:-1], hist_values, where="mid", label=component_labels.get(comp, comp),
                            linewidth=2, alpha=0.7, color=colors[i % len(colors)])
        handles.append(line)
    
    plt.xlabel(r"$m_{miss}^2$ (GeV$^2$)", fontsize=44)
#    plt.ylabel("Counts", fontsize=14)
    plt.ylabel("Normalised Counts", fontsize=34)
    plt.title(f"Bin {bin_index+1}")
    plt.legend(handles=handles, loc="best", fontsize=15, frameon=False)
    #save_path = os.path.join(output_dir, f"mmiss2_templates_bin{bin_index}.png")
    plt.tight_layout()
    save_path = os.path.join(output_dir, f"Normalised_mmiss2_templates_bin{bin_index}.png")
    print(f"Saving: {save_path}")
    plt.savefig(save_path, dpi=300)#, bbox_inches="tight")
    plt.close()

    ## Save each components sum to .json file
    os.makedirs('yields', exist_ok=True)
    output_json = os.path.join('yields', f"yields_bin{bin_index}.json")

    with open(output_json, 'w') as f:
        json.dump(yields, f, indent=4)

           


