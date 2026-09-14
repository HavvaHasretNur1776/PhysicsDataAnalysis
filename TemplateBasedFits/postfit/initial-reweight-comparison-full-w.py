#!/usr/bin/env python3
import ROOT
import json
import glob
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep  # LHCb style
import os
import sys


# Apply LHCb style
hep.style.use("LHCb2")

#ROOT.EnableImplicitMT()

# --------------------------------------------------
# Utilities
# --------------------------------------------------

# convert trigger weights to trigger bool decisions
'''def convert_trigger_weights_to_decisions(rdf):
    # gRandom->Uniform(0, 1) generates a random number between 0 and 1 for each event
    rdf = rdf.Define("rndm_val", "gRandom->Uniform(0, 1)")
    if "b0_L0Global_TIS" not in rdf.GetColumnNames():

        # Evaluate the condition entirely inside the RDataFrame graph
        rdf = rdf.Define("b0_L0Global_TIS", "b0_l0_global_tis_emu > rndm_val")
        rdf = rdf.Define("d0_L0HadronDecision_TOS", "d0_l0_hadron_tos_emu > rndm_val")
        #rdf = rdf.Define("d0_Hlt1TrackMVADecision_TOS", "d0_hlt1_trackmva_tos_emu")
        
    return rdf'''
def convert_trigger_weights_to_decisions(rdf):
    
    if "b0_L0Global_TIS" not in rdf.GetColumnNames():
        rdf = rdf.Define("rndm_val", "gRandom->Uniform(0, 1)")
        # Explicit double cast for NumPy compatibility
        rdf = rdf.Define("b0_L0Global_TIS", "(int)(b0_l0_global_tis_emu > rndm_val)")
        
    #if "d0_L0HadronDecision_TOS" not in rdf.GetColumnNames():
        rdf = rdf.Define("d0_L0HadronDecision_TOS", "(int)(d0_l0_hadron_tos_emu > rndm_val)")
        
    return rdf

def make1DHistos(data, bins, weights, binrange):
    if weights is None:
        hist, b = np.histogram(data, bins=bins, range=binrange)
        var = np.where(hist > 0, hist, 1.0)
    else:
        hist, b = np.histogram(data, bins=bins, range=binrange, weights=weights)
        var, _ = np.histogram(
            data, bins=bins, range=binrange, weights=weights**2
        )
        #avg_sqweight = np.mean(weights) ** 2
        avg_sqweight = np.mean(weights ** 2)
        var = np.where(hist > 0, var, avg_sqweight)

    return hist, np.sqrt(var), b


def makeDataHist(infiles, tree, branch, nbins, binrange):
    h_tot = None
    herr2_tot = None

    for infile in infiles:
        rdf = ROOT.RDataFrame(tree, infile)

        rdf = rdf.Filter("m2miss_reco_approx > -5 & m2miss_reco_approx < 6")

        # add eta
        if "ETA" in branch:
            part = branch.split("_")[0]
            pz_branch = f"{part}_PZ"
            pt_branch = f"{part}_PT"
            rdf = rdf.Define(f"{part}_ETA",f"asinh({pz_branch}/{pt_branch})")

        data = rdf.AsNumpy([branch])[branch]

        h, b = np.histogram(data, bins=nbins, range=binrange)
        err2 = np.where(h > 0, h, 1.0)

        if h_tot is None:
            h_tot = h
            herr2_tot = err2
        else:
            h_tot += h
            herr2_tot += err2

    return h_tot, np.sqrt(herr2_tot), b


# --------------------------------------------------
# Post-fit cocktail
# --------------------------------------------------

class PostFitCocktail:
    def __init__(self):
        self.samples = {}

    def addSample(self, name, infiles, tree="DecayTree", weightbranch=None):
        self.samples[name] = {
            "infiles": infiles,
            "tree": tree,
            "w": weightbranch,
        }

    def makePostFitHist(self, branch, nbins, binrange, fit_results):
        h = None
        herr2 = None

        for name, sample in self.samples.items():
            '''if name not in fit_results:
                print(f"[PostFitCocktail] Skipping {name}: not in fit_results")
                continue'''

            N = fit_results[name]["N"]
            Nerr = fit_results[name]["Nerr"]

            ih_tot = None
            iherr2_tot = None

            for infile in sample["infiles"]:
                rdf = ROOT.RDataFrame(sample["tree"], infile)
                rdf = rdf.Filter("m2miss_reco_approx > -5 & m2miss_reco_approx < 6")
                rdf = convert_trigger_weights_to_decisions(rdf)
                part = branch.split("_")[0]
                if "ETA" in branch and f"{part}_ETA" not in rdf.GetColumnNames():
                    #part = branch.split("_")[0]
                    pz_branch = f"{part}_PZ"
                    pt_branch = f"{part}_PT"
                    rdf = rdf.Define(f"{part}_ETA",f"asinh({pz_branch}/{pt_branch})")

                data = rdf.AsNumpy([branch])[branch]

                weights = None
                if sample["w"]:
                    #print(f"Using weights for sample {name}: {sample['w']}")
                    rdf = rdf.Define("total_weight", sample["w"])
                    weights = rdf.AsNumpy(["total_weight"])["total_weight"]
                    #weights = rdf.AsNumpy([sample["w"]])[sample["w"]]

                ih, iherr, b = make1DHistos(data, nbins, weights, binrange)

                if ih_tot is None:
                    ih_tot = ih
                    iherr2_tot = iherr**2
                else:
                    ih_tot += ih
                    iherr2_tot += iherr**2

            if ih_tot is None:
                continue

            norm = np.sum(ih_tot)
            #print(name, "raw sum:", np.sum(ih_tot))

            if norm == 0:
                continue

            ih_norm = ih_tot / norm
            ih_norm_err = np.sqrt(iherr2_tot) / norm

            ih_pf = N * ih_norm
            #ih_pf_var = (N * ih_norm_err) ** 2 + (ih_norm * Nerr) ** 2
            # Only use the statistical uncertainty of the MC templates for the bin errors
            ih_pf_var = (N * ih_norm_err) ** 2

            if h is None:
                h = ih_pf
                herr2 = ih_pf_var
            else:
                h += ih_pf
                herr2 += ih_pf_var
        
        return h, np.sqrt(herr2), b



# --------------------------------------------------
# Fit results loader
# --------------------------------------------------

def load_fullw_fit_results(json_pattern):
    fit_sum = {}

    for jf in sorted(glob.glob(json_pattern)):
        with open(jf) as f:
            data = json.load(f)

        for name, y in data["yields"].items():
            val = float(y["val"])
            err = float(y["err"])

            if name not in fit_sum:
                fit_sum[name] = {"N": 0.0, "Nerr2": 0.0}

            fit_sum[name]["N"] += val
            fit_sum[name]["Nerr2"] += err**2

    for name in fit_sum:
        fit_sum[name]["Nerr"] = np.sqrt(fit_sum[name]["Nerr2"])
        del fit_sum[name]["Nerr2"]

    return fit_sum


# --------------------------------------------------
# Configuration
# --------------------------------------------------

#EOS = "root://eoslhcb.cern.ch"
years = [2016, 2017, 2018]

file_map = {
    "nDstMu": ["b0_dstmunu_up", "b0_dstmunu_down"],
    "nDstTau": ["b0_dsttaunu_up", "b0_dsttaunu_down"],

    "nDstDX": [
        "b+_dstdx_up","b+_dstdx_down",
        "b0_dstdx_up","b0_dstdx_down"
    ],

    "nDststHigher": ["b0_dststmu_heavy_up", "b+_dststmu_heavy_up", "b0_dststmu_heavy_down", "b+_dststmu_heavy_down"],

    "nb_dstst": [
        "b0_dststmu_up", "b0_dststmu_down",
        "b+_dststmu_up", "b+_dststmu_down"
    ],

    "nbs_dstst": [

        "bs_dststmu_up", "bs_dststmu_down"
    ],
    "nWSPi": ["data_WSPi_up", "data_WSPi_down"],
    "nWSMu": ["data_WSMu_up", "data_WSMu_down"]
}
base_path_mc = "/eos/lhcb/user/d/dvicoben/bd2dstmunu_run2/Data/MCTracker"
base_path_data = "/eos/lhcb/user/d/dvicoben/bd2dstmunu_run2/Data/Collision"

cocktail = PostFitCocktail()
for comp, subpaths in file_map.items():
    #print(f"Adding component: {comp} with subpaths: {subpaths}")
    if comp in ["nWSPi", "nWSMu"]: # data-driven samples do not use weights
        infiles = [
            f"{base_path_data}/{y}/selected/{sp}_selected.root"
            for y in years
            for sp in subpaths
        ]
        cocktail.addSample(comp, infiles)

    # Add weights: PID+L0+tracking+SPD+bkin+hammer weights for MC samples, with different hammer weights for BGL vs BLR schemes
    elif comp in ["nDstMu", "nDstTau"]:
        infiles = [
            f"{base_path_mc}/{y}/selected/{sp}_selected.root"
            for y in years
            for sp in subpaths
        ]
        total_weight_BGL = "PIDCalibEff*w_trk_mu*w_trk_k*w_trk_pi*w_trk_spi*w_l0_tosORtis_emu*w_spdcorr*w_bkin_tot*SchemeBGL"
        cocktail.addSample(comp, infiles, weightbranch=total_weight_BGL)#
    
    elif comp in ["nb_dstst", "nbs_dstst"]:
        infiles = [
            f"{base_path_mc}/{y}/selected/{sp}_selected.root"
            for y in years
            for sp in subpaths
        ]
        total_weight_BLR = "PIDCalibEff*w_trk_mu*w_trk_k*w_trk_pi*w_trk_spi*w_l0_tosORtis_emu*w_spdcorr*w_bkin_tot*SchemeBLR"
        cocktail.addSample(comp, infiles, weightbranch=total_weight_BLR)#
    
    else: # non-hammer MC samples
        infiles = [
            f"{base_path_mc}/{y}/selected/{sp}_selected.root"
            for y in years
            for sp in subpaths

        ]
        total_weight_non_hammer = "PIDCalibEff*w_trk_mu*w_trk_k*w_trk_pi*w_trk_spi*w_l0_tosORtis_emu*w_spdcorr*w_bkin_tot"

        cocktail.addSample(comp, infiles, weightbranch=total_weight_non_hammer)#
      
    print(f"Sample added to cocktail: {comp} ")
# --------- End of cocktail setup ----------------------------------


# Set up a cocktail without b-kin weights for comparison
cocktail_noWeight = PostFitCocktail()
for comp, subpaths in file_map.items():
    if comp in ["nWSPi", "nWSMu"]: # data-driven samples do not use track weights
        infiles = [
            f"{base_path_data}/{y}/selected/{sp}_selected.root"
            for y in years
            for sp in subpaths
        ]
        cocktail_noWeight.addSample(comp, infiles)

    elif comp in ["nDstMu", "nDstTau"]:
        infiles = [
            f"{base_path_mc}/{y}/selected/{sp}_selected.root"
            for y in years
            for sp in subpaths
        ]
        total_weight_BGL_noBkin = "PIDCalibEff*w_trk_mu*w_trk_k*w_trk_pi*w_trk_spi*w_l0_tosORtis_emu*w_spdcorr*SchemeBGL"
        cocktail_noWeight.addSample(comp, infiles, weightbranch=total_weight_BGL_noBkin)#
    
    elif comp in ["nb_dstst", "nbs_dstst"]:
        infiles = [
            f"{base_path_mc}/{y}/selected/{sp}_selected.root"
            for y in years
            for sp in subpaths
        ]
        total_weight_BLR_noBkin = "PIDCalibEff*w_trk_mu*w_trk_k*w_trk_pi*w_trk_spi*w_l0_tosORtis_emu*w_spdcorr*SchemeBLR"
        cocktail_noWeight.addSample(comp, infiles, weightbranch=total_weight_BLR_noBkin)#
    
    else: # non-hammer MC samples
        infiles = [
            f"{base_path_mc}/{y}/selected/{sp}_selected.root"
            for y in years
            for sp in subpaths
        ]
        total_weight_non_hammer_noBkin = "PIDCalibEff*w_trk_mu*w_trk_k*w_trk_pi*w_trk_spi*w_l0_tosORtis_emu*w_spdcorr"

        cocktail_noWeight.addSample(comp, infiles, weightbranch=total_weight_non_hammer_noBkin)#
      
    print(f"Sample added to no-weight cocktail: {comp} ")

data_files = [
    f"{base_path_data}/{y}/selected/data_RS_{m}_selected.root"
    for y in years
    for m in ["up", "down"]
]

fit_results = load_fullw_fit_results(
    f"../fullw_fitted_yields_v2/fittedYields_fullq2.json"
)


# --------------------------------------------------
# Plot 
# --------------------------------------------------
# Create output directory for plots
output_dir = "comparison-initial-reweighting"


os.makedirs(output_dir, exist_ok=True)

# Variables to plot
branch_info = { 
    # B0
    "b0_PT": (100, (0.0, 40000.0), r"$p_{T}(B^{0})\,[\mathrm{MeV}/c]$"),
    #"b0_P": (100, (10000.0, 200000.0), r"$p(B^{0})\,[\mathrm{MeV}/c]$"),
    "b0_ETA": (24, (2.0, 6.0), r"$\eta(B^{0})$"),
    "b0_IP_OWNPV": (50, (0.0, 2.0), r"$\mathrm{IP}(B^{0})\,[\mathrm{mm}]$"),
    "b0_FDCHI2_OWNPV": (50, (0.0, 10000.0), r"$\mathrm{FD}\chi^{2}(B^{0})$"),
    "b0_FD_OWNPV": (50, (0.0, 50.0), r"$\mathrm{FD}(B^{0})\,[\mathrm{mm}]$"),
    "b0_OWNPV_NDOF": (75, (0.0, 300.0), "PV/ndf"),
    "nTracks": (20, (0.0, 400.0), "nTracks"),
    "b0_L0Global_TIS": (2, (-0.5, 1.5), "L0Global_TIS"),
    "d0_L0HadronDecision_TOS": (2, (-0.5, 1.5), "L0Hadron_TOS"),

    # Muon
    "mu_PT": (80, (1000.0, 20000.0), r"$p_{T}(\mu)\,[\mathrm{MeV}/c]$"),
    #"mu_P": (100, (10000.0, 100000.0), r"$p(\mu)\,[\mathrm{MeV}/c]$"),
    #"mu_ETA": (24, (2.0, 6.0), r"$\eta(\mu)$"),

    # D0
    "d0_PT": (100, (2000.0, 15000.0), r"$p_{T}(D^{0})\,[\mathrm{MeV}/c]$"),
    #"d0_P": (80, (10000.0, 150000.0), r"$p(D^{0})\,[\mathrm{MeV}/c]$"),
    #"d0_ETA": (24, (2.0, 6.0), r"$\eta(D^{0})$"),
    #"d0_IP_OWNPV": (50, (0.0, 2.0), r"$\mathrm{IP}(D^{0})\,[\mathrm{mm}]$"),
    #"d0_FDCHI2_OWNPV": (50, (0.0, 100000.0), r"$\mathrm{FD}\chi^{2}(D^{0})$"),

    # D*
    #"dst_FDCHI2_OWNPV": (50, (0.0, 1000.0), r"$\mathrm{FD}\chi^{2}(D^{*})$"),

    # Kaon
    "k_PT": (80, (900.0, 10000.0), r"$p_{T}(K)\,[\mathrm{MeV}/c]$"),
    #"k_P": (100, (0.0, 100000.0), r"$p(K)\,[\mathrm{MeV}/c]$"),
    #"k_ETA": (24, (2.0, 6.0), r"$\eta(K)$"),
    #"k_IP_OWNPV": (50, (0.0, 5.0), r"$\mathrm{IP}(K)\,[\mathrm{mm}]$"),

    # Pion
    "pi_PT": (80, (900.0, 9000.0), r"$p_{T}(\pi)\,[\mathrm{MeV}/c]$"),
    #"pi_P": (100, (0.0, 100000.0), r"$p(\pi)\,[\mathrm{MeV}/c]$"),
    #"pi_ETA": (24, (2.0, 6.0), r"$\eta(\pi)$"),
    #"pi_IP_OWNPV": (50, (0.0, 5.0), r"$\mathrm{IP}(\pi)\,[\mathrm{mm}]$")

}


# Loop over branches to plot
for branch, (nbins, binrange, xlabel) in branch_info.items():

    # Create figure per branch
    fig, (ax, axr) = plt.subplots(
        2, 1,
        figsize=(10, 9),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True
    )

    # Post-fit MC
    h_mc, herr_mc, b = cocktail.makePostFitHist(
        branch, nbins, binrange, fit_results
    )

    # Post-fit MC without b-kin weights for comparison
    h_mc_noWeight, herr_mc_noWeight, b = cocktail_noWeight.makePostFitHist(
        branch, nbins, binrange, fit_results
    )

    # Data
    h_data, herr_data, _ = makeDataHist(
        data_files, "DecayTree", branch, nbins, binrange
    )

    bc = 0.5 * (b[1:] + b[:-1])

    # --- Main plot ---
    
    # normalize everything to unit area for shape comparison
    h_mc_norm = h_mc / np.sum(h_mc) if np.sum(h_mc) > 0 else h_mc
    herr_mc_norm = herr_mc / np.sum(h_mc) if np.sum(h_mc) > 0 else herr_mc

    h_mc_noWeight_norm = h_mc_noWeight / np.sum(h_mc_noWeight) if np.sum(h_mc_noWeight) > 0 else h_mc_noWeight
    herr_mc_noWeight_norm = herr_mc_noWeight / np.sum(h_mc_noWeight) if np.sum(h_mc_noWeight) > 0 else herr_mc_noWeight

    h_data_norm = h_data / np.sum(h_data) if np.sum(h_data) > 0 else h_data
    herr_data_norm = herr_data / np.sum(h_data) if np.sum(h_data) > 0 else herr_data

    #ax.stairs(h_mc_norm, b, color="blue")
    ax.errorbar(#post-fit
        #bc, h_mc, yerr=herr_mc,
        bc, h_mc_norm, yerr=herr_mc_norm,
        fmt="o", markersize=8, capsize=4, color="blue", label="After reweighting"
    )


    ### no b-kin weights for comparison
    #ax.stairs(h_mc_noWeight_norm, b, color="green", linestyle="--")
    ax.errorbar(#post-fit no b-kin weights
        #bc, h_mc_noWeight, yerr=herr_mc_noWeight,
        bc, h_mc_noWeight_norm, yerr=herr_mc_noWeight_norm,
        fmt="o", markersize=8, capsize=4, color="green", label="Before reweighting"
    )

    ax.errorbar(#data
        #bc, h_data, yerr=herr_data,
        bc, h_data_norm, yerr=herr_data_norm,
        fmt="o", color="black", markersize=8,
        capsize=4, label="Data"
    )


    ax.set_ylabel("Events")
    ax.set_title(f"$1.\leq w < 1.52$")
    ax.legend(fontsize=17)

    # --- Ratio plot ---
    ratio = np.zeros_like(h_data, dtype=float)
    ratio_err = np.zeros_like(h_data, dtype=float)

    mask = (h_mc > 0) & (h_data > 0)
    ratio[mask] = h_data[mask] / h_mc[mask]
    ratio_err[mask] = ratio[mask] * np.sqrt(
        (herr_data[mask] / h_data[mask])**2 +
        (herr_mc[mask] / h_mc[mask])**2
    )

    axr.axhline(1.0, color="gray", linestyle="--")
    axr.errorbar(
        bc, ratio, yerr=ratio_err,
        fmt="o", color="blue", markersize=8,
        capsize=4
    )

    #plot also the ratio for no b-kin weights for comparison
    ratio_noWeight = np.zeros_like(h_data, dtype=float)
    ratio_noWeight_err = np.zeros_like(h_data, dtype=float) 
    mask_noWeight = (h_mc_noWeight > 0) & (h_data > 0)
    ratio_noWeight[mask_noWeight] = h_data[mask_noWeight] / h_mc_noWeight[mask_noWeight]
    ratio_noWeight_err[mask_noWeight] = ratio_noWeight[mask_noWeight] * np.sqrt(
        (herr_data[mask_noWeight] / h_data[mask_noWeight])**2 +
        (herr_mc_noWeight[mask_noWeight] / h_mc_noWeight[mask_noWeight])**2
    )
    axr.errorbar(
        bc, ratio_noWeight, yerr=ratio_noWeight_err,
        fmt="o", color="green", markersize=8,
        capsize=4
    )

    axr.set_xlabel(xlabel, fontsize=56)
    axr.set_ylabel("Data / Model")
    axr.set_xlim(binrange)
    axr.set_ylim(0.5, 1.5)
#    axr.set_ylim(-2.5, 2.5)

    plt.tight_layout()
    output_plot = f"{output_dir}/{branch}.png"
    plt.savefig(output_plot, dpi=300)
    plt.close(fig)

    print(f"Plot saved as {output_plot}")


