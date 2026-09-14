# Run on lb-conda Semilep/hammer/2023-07-15 bash
import numpy as np
import uproot
import os
import ROOT
import sys
import matplotlib.pyplot as plt
import mplhep as hep  # LHCb style

# Apply LHCb style
hep.style.use("LHCb2")

#year = sys.argv[1]  # e.g., '2016', '2017', '2018'
#data_path = f"/eos/lhcb/user/h/havva/q2Measurement-BtoDstMuNu/MCsamples/TrackerOnlyNTuples_Maryland/Angular_RDst_muonic_Run2/{year}"
data_types = [
    "b0_dstmunu","b0_dsttau", 
    "b0_dststtau", "b+_dststtau",
    "b+_dstdx", "b0_dstdx",
    "b0_dststmu_415", "b0_dststmu_10413", "b0_dststmu_20413",
    "b+_dststmu_425", "b+_dststmu_10423", "b+_dststmu_20423",
    "bs_dststmu_10433", "bs_dststmu_435",
    "b+_dststmu_heavy", "b0_dststmu_heavy"
]

'''output_dir = f"weighted_histograms_pid+l0+hlt1_20bins/{year}/"
os.makedirs(output_dir, exist_ok=True)
output_dir_plot=f'weighted_histograms_pid+l0+hlt1_20bins/{year}/plot_{year}'
os.makedirs(output_dir_plot, exist_ok=True)
'''
# Define LaTeX-style labels
component_labels = {
    "b0_dstmunu": r"$B^0 \to D^{*+} \mu^- \nu_\mu$",
    "b0_dsttau": r"$B^0 \to D^{*+} \tau^- \nu_\tau$",
    "b0_dststtau": r"$B^0 \to D^{**+} \tau^- \nu_\tau$",
    "b+_dststtau": r"$B^+ \to D^{**+} \tau^- \nu_\tau$",
    "b+_dstdx": r"$B^+ \to D^{*+}D(X_c \to \mu \nu X)$",
    "b0_dstdx": r"$B \to D^{*}D(X_c \to \mu \nu X)$",
    "b0_dststmu_415": r"$B^0 \to D_2^*(2460)^+ \mu^- \nu_\mu$",
    "b+_dststmu_425": r"$B^+ \to D_2^*(2460)^0 \mu^- \nu_\mu$",
    "b0_dststmu_10413": r"$B^0 \to D_1(2420)^+ \mu^- \nu_\mu$",
    "b0_dststmu_20413": r"$B^0 \to D_1(H)^+ \mu^- \nu_\mu$",
    "b+_dststmu_10423": r"$B^+ \to D_1(2420)^0 \mu^- \nu_\mu$",
    "b+_dststmu_20423": r"$B^+ \to D_1(2430)^0 \mu^- \nu_\mu$",
    "bs_dststmu_435": r"$B_s \to D_{s2}^*(2573)^+ \mu^- \nu_\mu$",
    "bs_dststmu_10433": r"$B_s \to D_{s1}(2536)^+ \mu^- \nu_\mu$",
    "b0_dststmu_heavy": r"$B^0 \to D^{**}(Higher) \mu^- \nu_\mu$",
    "b+_dststmu_heavy": r"$B^+ \to D^{**}(Higher) \mu^- \nu_\mu$"

}
# Define colors
colors = plt.cm.tab20.colors

# Define binning masks and bin edges
def get_masks(w_Pred):
    # 1 <= w <= 1.5
    return[
        w_Pred <= 1.050,
        (w_Pred > 1.050) & (w_Pred <= 1.100),
        (w_Pred > 1.100) & (w_Pred <= 1.150),
        (w_Pred > 1.150) & (w_Pred <= 1.200),
        (w_Pred > 1.200) & (w_Pred <= 1.250),
        (w_Pred > 1.250) & (w_Pred <= 1.300),
        (w_Pred > 1.300) & (w_Pred <= 1.350),
        (w_Pred > 1.350) & (w_Pred <= 1.400),
        (w_Pred > 1.400) & (w_Pred <= 1.450),
        w_Pred > 1.450 
    ]
    # 10 bins: 1.0, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, >1.45
'''custom_mmiss2_bins = [
    [-4e6, 5.5e6],    # Bin 1
    [-5e6, 6e6],      # Bin 2
    [-3e6, 5.5e6],    # Bin 3
    [-3e6, 5e6],      # Bin 4
    [-3e6, 5e6],      # Bin 5
    [-3e6, 5e6],      # Bin 6
    [-2.8e6, 5e6],    # Bin 7
    [-2e6, 5e6],      # Bin 8
    [-2e6, 4.1e6],    # Bin 9
    [-2e6, 3e6]       # Bin 10
]'''
'''mmiss2_bin_edges = [
    np.linspace(-4., 5.5, 21),    # Bin 1
    np.linspace(-5., 6., 21),   # Bin 2
    np.linspace(-3., 5.5, 21),      # Bin 3
    np.linspace(-3., 5.5, 21),      # Bin 4
    np.linspace(-3., 5., 21),    # Bin 5
    np.linspace(-3., 5., 21),    # Bin 6
    np.linspace(-2.8, 5., 21), # Bin 7
    np.linspace(-2., 5., 21),  # Bin 8
    np.linspace(-1.0, 4.1, 21),  # Bin 9
    np.linspace(-2., 3., 21),  # Bin 10
]'''
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
#bin_edges = np.linspace(-5e6, 6e6, 21)
bin_edges = np.linspace(-5., 6., 21)
#template for all w bins
for year in ['2016','2017','2018']:
    data_path = f"./{year}"
    output_dir = f"pid_weighted_histograms/{year}/"
    os.makedirs(output_dir, exist_ok=True)

    for component in data_types:
#        path = f'{data_path}/{component}_BestVariables.root'
        path = f'{data_path}/{component}_BestVariables.root'
        file_data = ROOT.TFile.Open(path)
        tree_data = file_data.Get("DecayTree")
        df = ROOT.RDataFrame(tree_data)

        # Extract variables - this type of extraction avoids memory crashes-
        mmiss2 = np.array(df.AsNumpy(columns=["mmiss2rec_brf"])["mmiss2rec_brf"])
        w_Pred = np.array(df.AsNumpy(columns=["w_Pred"])["w_Pred"])
        pid_eff = np.array(df.AsNumpy(columns=["pid_eff"])["pid_eff"])
        weights = pid_eff
        counts, _ = np.histogram(mmiss2, bins=bin_edges, weights=weights)

        # Save overall histogram
        with open(f"{output_dir}/{component}.dat", "w") as f:
            for c in counts:
                f.write(f"{c}\n")
        print(f"{component} - Overall histogram processed successfully.")

##templates per w bin
for year in ['2016','2017','2018']:
    data_path = f"./{year}"
    output_dir = f"pid_weighted_histograms/{year}/"
    output_dir_plot=f'pid_weighted_histograms/{year}/plot_{year}'
    os.makedirs(output_dir_plot, exist_ok=True)
    for i, mmiss2_range in enumerate(custom_mmiss2_bins, 1):
        handles = []#for plot legend
        plt.figure(figsize=(8, 6))

        for j,component in enumerate(data_types):#for coloring correctly
            path = f'{data_path}/{component}_BestVariables.root'
            file_data = ROOT.TFile.Open(path)
            tree_data = file_data.Get("DecayTree")
            df = ROOT.RDataFrame(tree_data)

            # Extract variables - this type of extraction avoids memory crashes-
            mmiss2 = np.array(df.AsNumpy(columns=["mmiss2rec_brf"])["mmiss2rec_brf"])
            w_Pred = np.array(df.AsNumpy(columns=["w_Pred"])["w_Pred"])
            pid_eff = np.array(df.AsNumpy(columns=["pid_eff"])["pid_eff"])
            weights = pid_eff  
            # apply bin mask
            mask = get_masks(w_Pred)[i-1]   # since get_masks gives list of masks
            mmiss2_bin = mmiss2[mask]
            weights_bin = weights[mask]
            bin_edges = np.linspace(mmiss2_range[0], mmiss2_range[1], 21)
            counts_bin, _ = np.histogram(mmiss2_bin, bins=bin_edges, weights=weights_bin)
            # Save per-bin histogram
            with open(f"{output_dir}/{component}_bin{i}.dat", "w") as f:
                for c in counts_bin:
                    f.write(f"{c}\n")
            print(f"{component} - Bin {i} processed successfully.")

            ##Plot each comp in each w bin
            yield_val = np.sum(counts_bin)
            label_with_yield = f"{component_labels.get(component, component)} ({yield_val:.1f})"  # Append yield
            counts_bin = counts_bin / yield_val  # Normalize counts
            line, = plt.step(bin_edges[:-1], counts_bin, linewidth=3, alpha=0.7, label=component_labels.get(component, component),color=colors[j % len(colors)])
            #line, = plt.step(bin_edges[:-1], counts_bin, linewidth=2, alpha=0.7, label=label_with_yield,color=colors[j % len(colors)])
            handles.append(line)
            plt.legend(handles=handles, loc="upper left", fontsize=10, frameon=True)

        plt.xlabel(r"$m_{miss}^2$ (GeV$^2$)", fontsize=14)
        plt.ylabel("Normalised Counts", fontsize=14)
        plt.title(f"Bin_{i}")
        plot_path = os.path.join(output_dir_plot, f"mmiss2_all_Normalised_components_bin{i}.png")
        xlow, xhigh = bin_edges[0], bin_edges[-1]
        plt.xlim(xlow, xhigh)
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Saved: {plot_path}")







