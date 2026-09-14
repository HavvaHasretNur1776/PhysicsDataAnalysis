#!/usr/bin/env python3
import ROOT
import json
import glob
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep  # LHCb style
import os
import sys

# Run ROOT in batch mode to prevent canvas pop-ups
ROOT.gROOT.SetBatch(True)

# Apply LHCb style
hep.style.use("LHCb2")

# --------------------------------------------------
# Utilities
# --------------------------------------------------

# convert trigger weights to trigger bool decisions
def convert_trigger_weights_to_decisions(rdf):
    # gRandom->Uniform(0, 1) generates a random number between 0 and 1 for each event
    rdf = rdf.Define("rndm_val", "gRandom->Uniform(0, 1)")
    if "b0_L0Global_TIS" not in rdf.GetColumnNames():

        # Evaluate the condition entirely inside the RDataFrame graph
        rdf = rdf.Define("b0_L0Global_TIS", "b0_l0_global_tis_emu > rndm_val")
        rdf = rdf.Define("d0_L0HadronDecision_TOS", "d0_l0_hadron_tos_emu > rndm_val")
        #rdf = rdf.Define("d0_Hlt1TrackMVADecision_TOS", "d0_hlt1_trackmva_tos_emu")
        
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
        rdf.Filter("m2miss_reco_approx > -5.0 && m2miss_reco_approx < 6.0")
        # need to match with old dltaM cut since these samples chaged!!
        ## NEW deltaM= 140-154  in samples
        cutdst_M = '(dst_M - d0_M) < 147 and (dst_M - d0_M) > 144' 
        rdf = rdf.Filter(cutdst_M)
        # add eta
        if "ETA" in branch:
            part = branch.split("_")[0]
            pz_branch = f"{part}_PZ"
            pt_branch = f"{part}_PT"
            if f"{part}_ETA" not in rdf.GetColumnNames():
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
            N = fit_results[name]["N"]
            Nerr = fit_results[name]["Nerr"]

            ih_tot = None
            iherr2_tot = None

            for infile in sample["infiles"]:
                rdf = ROOT.RDataFrame(sample["tree"], infile)
                rdf.Filter("m2miss_reco_approx > -5.0 && m2miss_reco_approx < 6.0")
                # apply the same deltaM cut as in makeDataHist
                cutdst_M = '(dst_M - d0_M) < 147 and (dst_M - d0_M) > 144'
                rdf = rdf.Filter(cutdst_M)
                rdf = convert_trigger_weights_to_decisions(rdf)

                # --- NEW ROBUST ETA DEFINITION ---
                # We need to check both the branch name AND the weight string
                # for any particle_ETA variables.
                weight_str = sample["w"] if sample["w"] else ""
                
                # Particles that might have an ETA defined in your weights
                potential_particles = ["b0", "d0", "mu", "pi", "k", "spi"]

                for part in potential_particles:
                    eta_var = f"{part}_ETA"
                    # If ETA is needed for the plot OR mentioned in the weights
                    if (eta_var in branch) or (weight_str and eta_var in weight_str):
                        if eta_var not in rdf.GetColumnNames():
                            pz = f"{part}_PZ"
                            pt = f"{part}_PT"
                            # Extra safety check: ensure PZ/PT exist before defining ETA
                            if pz in rdf.GetColumnNames() and pt in rdf.GetColumnNames():
                                rdf = rdf.Define(eta_var, f"asinh({pz}/{pt})")
                
                # Now that all columns (ETA, trigger, etc.) are defined, Define the weight
                if sample["w"]:
                    rdf = rdf.Define("total_weight", sample["w"])
                    data_dict = rdf.AsNumpy([branch, "total_weight"])
                    data = data_dict[branch]
                    weights = data_dict["total_weight"]
                else:
                    data = rdf.AsNumpy([branch])[branch]
                    weights = None

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

            if norm == 0:
                continue

            ih_norm = ih_tot / norm
            ih_norm_err = np.sqrt(iherr2_tot) / norm

            ih_pf = N * ih_norm
#            ih_pf_var = (N * ih_norm_err) ** 2 + (ih_norm * Nerr) ** 2
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

years = [2016, 2017, 2018]

file_map = {
    "nDstMu": ["b0_dstmunu_up", "b0_dstmunu_down"],
    "nDstTau": ["b0_dsttaunu_up", "b0_dsttaunu_down"],
    "nDstDX": ["b+_dstdx_up","b+_dstdx_down", "b0_dstdx_up","b0_dstdx_down"],
    "nDststHigher": ["b0_dststmu_heavy_up", "b+_dststmu_heavy_up", "b0_dststmu_heavy_down", "b+_dststmu_heavy_down"],
    "nb_dstst": ["b0_dststmu_up", "b0_dststmu_down", "b+_dststmu_up", "b+_dststmu_down"],
    "nbs_dstst": ["bs_dststmu_up", "bs_dststmu_down"],
    "nWSPi": ["data_WSPi_up", "data_WSPi_down"],
    "nWSMu": ["data_WSMu_up", "data_WSMu_down"]
}

base_path_mc = "/eos/lhcb/user/d/dvicoben/bd2dstmunu_run2/Data/MCTracker"
base_path_data = "/eos/lhcb/user/d/dvicoben/bd2dstmunu_run2/Data/Collision"

cocktail = PostFitCocktail()
for comp, subpaths in file_map.items():
    if comp in ["nWSPi", "nWSMu"]:
        infiles = [f"{base_path_data}/{y}/selected/{sp}_selected.root" for y in years for sp in subpaths]
        cocktail.addSample(comp, infiles)
    elif comp in ["nDstMu", "nDstTau"]:
        infiles = [f"{base_path_mc}/{y}/selected/{sp}_selected.root" for y in years for sp in subpaths]
        total_weight_BGL = "PIDCalibEff*w_trk_mu*w_trk_k*w_trk_pi*w_trk_spi*w_l0_tosORtis_emu*w_spdcorr*w_bkin_tot*SchemeBGL"
        cocktail.addSample(comp, infiles, weightbranch=total_weight_BGL)
    elif comp in ["nb_dstst", "nbs_dstst"]:
        infiles = [f"{base_path_mc}/{y}/selected/{sp}_selected.root" for y in years for sp in subpaths]
        total_weight_BLR = "PIDCalibEff*w_trk_mu*w_trk_k*w_trk_pi*w_trk_spi*w_l0_tosORtis_emu*w_spdcorr*w_bkin_tot*SchemeBLR"
        cocktail.addSample(comp, infiles, weightbranch=total_weight_BLR)
    else: 
        infiles = [f"{base_path_mc}/{y}/selected/{sp}_selected.root" for y in years for sp in subpaths]
        total_weight_non_hammer = "PIDCalibEff*w_trk_mu*w_trk_k*w_trk_pi*w_trk_spi*w_l0_tosORtis_emu*w_spdcorr*w_bkin_tot"
        cocktail.addSample(comp, infiles, weightbranch=total_weight_non_hammer)
      
    print(f"Sample added to cocktail: {comp} ")

data_files = [f"{base_path_data}/{y}/selected/data_RS_{m}_selected.root" for y in years for m in ["up", "down"]]

fit_results = load_fullw_fit_results(
    f"../fullw_fitted_yields_v2/fittedYields_fullq2.json"
)

# --------------------------------------------------
# Sequential 3D Reweighting Engine
# --------------------------------------------------

# Example configuration for reweighting stages. Adjust variables, bins, and ranges as needed.
reweighting_stages = [
    
    {
        "vars": ["nTracks", "b0_OWNPV_NDOF", "b0_L0Global_TIS"],
        "bins": [10, 10, 2], 
        "ranges": [(0.0, 400), (0.0, 200.0), (0.0, 2.0)],
        "file_name": "nTracks_b0_OWNPV_NDOF_b0_L0Global_TIS"

    },
   
    # makes d0-pt and d0/k/pi-IP worse comparing to the overall stages
    #{
    #    "vars":["d0_PT","log(d0_IPCHI2_OWNPV)", "d0_IPCHI2_OWNPV/d0_IP_OWNPV"],
    #    "bins":[20, 20, 20],
    #    "ranges":[(2000.0, 15000.0), (2.0, 10.0), (5.0, 110.0)],
    #    "file_name": "d0_PT_log_d0_IPCHI2_OWNPV_d0_IPCHI2_over_d0_IP_OWNPV"

    #},
    
    {
        "vars":["d0_PT","d0_ETA", "b0_L0Global_TIS"],
        "bins":[20, 6, 2],
        "ranges":[(1000.0, 20000.0), (0.0, 5.0), (0.0, 2.0)],
        "file_name": "d0_PT_d0_ETA_b0_L0Global_TIS"
    },
    {
        "vars": ["k_PT", "k_ETA", "b0_L0Global_TIS"],
        "bins": [10, 6, 2],
        "ranges": [(1000.0, 8000.0), (0.0, 5.0), (0.0, 2.0)],
        "file_name": "k_PT_k_ETA_b0_L0Global_TIS"
    },

    {
        "vars": ["pi_PT", "pi_ETA", "b0_L0Global_TIS"],
        "bins": [10, 6, 2],     
        "ranges": [(1000.0, 8000.0), (0.0, 5.0), (0.0, 2.0)],
        "file_name": "pi_PT_pi_ETA_b0_L0Global_TIS"
    },
    {
        "vars": ["mu_PT", "mu_ETA", "b0_L0Global_TIS"],
        "bins": [20, 6, 2],
        "ranges": [(1000.0, 20000.0), (0.0, 5.0), (0.0, 2.0)],
        "file_name": "mu_PT_mu_ETA_b0_L0Global_TIS"
    },
    
    {
        "vars": ["log(1-d0_DIRA_OWNPV)", "b0_L0Global_TIS", "b0_L0Global_TIS"],
        "bins": [40, 2, 2],
        "ranges": [(-14.0, 8.4), (0.0, 2.0), (0.0, 2.0)],
        "file_name": "log1m_d0_DIRA_OWNPV_b0_L0Global_TIS"

    }


   
]


def apply_sequential_reweighting(data_files, cocktail, fit_results, stages):
    # Track accumulated weight strings
    current_weights = {
        name: (samp["w"] if samp["w"] else "1.0") 
        for name, samp in cocktail.samples.items()
    }
    
    for i, stage in enumerate(stages):
        print(f"[INFO] Running 3D Reweighting Stage {i+1} using {stage['vars']}")
        vars_3d = stage["vars"]
        bins_3d = stage["bins"]
        ranges_3d = stage["ranges"]

        def prep_rdf(rdf, name=None, is_mc=False):
            
            if is_mc:
                rdf = convert_trigger_weights_to_decisions(rdf)

            # Get the weight string for this specific sample to check for dependencies
            # If name is None (like for Data), we just check the current stage vars
            weight_str = current_weights.get(name, "") if name else ""
            
            # 1. Define ALL potential ETAs if they appear in the weight OR the current stage
            potential_parts = ["b0", "d0", "mu", "pi", "k", "spi"]
            for part in potential_parts:
                eta_var = f"{part}_ETA"
                # Check if it's needed for this stage OR was used in a previous stage
                if any(eta_var in v for v in vars_3d) or (eta_var in weight_str):
                    if eta_var not in rdf.GetColumnNames():
                        rdf = rdf.Define(eta_var, f"asinh({part}_PZ/{part}_PT)")

            # 2. Now handle the current stage's special variables (logs, ratios, etc.)
            safe_vars = []
            for v in vars_3d:
                if "ETA" in v:
                    safe_vars.append(v) # Already defined above
                elif "log(1-" in v:
                    base_var = v.replace("log(1-", "").replace(")", "")
                    safe_name = f"log1m_{base_var}"
                    if safe_name not in rdf.GetColumnNames():
                        rdf = rdf.Define(safe_name, f"log(1 - {base_var})")
                    safe_vars.append(safe_name)
                elif "log(" in v:
                    base_var = v.replace("log(", "").replace(")", "")
                    safe_name = f"log_{base_var}"
                    if safe_name not in rdf.GetColumnNames():
                        rdf = rdf.Define(safe_name, f"log({base_var})")
                    safe_vars.append(safe_name)
                elif "/" in v:
                    num, denom = v.split("/")
                    safe_name = f"{num}_over_{denom}"
                    if safe_name not in rdf.GetColumnNames():
                        rdf = rdf.Define(safe_name, f"(sqrt({num}))/({denom})")
                    safe_vars.append(safe_name)
                else:
                    safe_vars.append(v)
                    
            return rdf, safe_vars
        
        
        # 1. Fill Data 3D Histogram
        rdf_data_full, safe_vars = prep_rdf(ROOT.RDataFrame("DecayTree", tuple(data_files)), is_mc=False)
        
        # Apply the specific region cut ONLY for filling the histogram
        rdf_data_sig = rdf_data_full.Filter("m2miss_reco_approx < 0.5")
        
        model_data = ROOT.RDF.TH3DModel(f"h_data_stage{i}", "Data", 
                                        bins_3d[0], ranges_3d[0][0], ranges_3d[0][1],
                                        bins_3d[1], ranges_3d[1][0], ranges_3d[1][1],
                                        bins_3d[2], ranges_3d[2][0], ranges_3d[2][1])

        h3_data_ptr = rdf_data_sig.Histo3D(model_data, safe_vars[0], safe_vars[1], safe_vars[2])

        # 2. Fill Post-Fit MC 3D Histograms
        # 2. Fill Post-Fit MC 3D Histograms
        mc_hist_ptrs = {}
        mc_norms_nocut = {} # Dictionary to hold the pre-cut integrals

        for name, sample in cocktail.samples.items():
            rdf_mc_full, _ = prep_rdf(ROOT.RDataFrame(sample["tree"], tuple(sample["infiles"])), name=name, is_mc=True)
            rdf_mc_full = rdf_mc_full.Define("current_w", current_weights[name])
            
            # Action 1: Get the total sum of weights in the FULL fit region
            mc_norms_nocut[name] = rdf_mc_full.Sum("current_w")

            # Action 2: Filter for the reweighting region and book the histogram
            rdf_mc_sig = rdf_mc_full.Filter("m2miss_reco_approx < 0.5")
            
            model_mc = ROOT.RDF.TH3DModel(f"h_{name}_stage{i}", name, 
                                          bins_3d[0], ranges_3d[0][0], ranges_3d[0][1],
                                          bins_3d[1], ranges_3d[1][0], ranges_3d[1][1],
                                          bins_3d[2], ranges_3d[2][0], ranges_3d[2][1])
            
            mc_hist_ptrs[name] = rdf_mc_sig.Histo3D(model_mc, safe_vars[0], safe_vars[1], safe_vars[2], "current_w")
        
        # Aggregate MC
        h3_data = h3_data_ptr.GetValue()
        h3_mc_total = h3_data.Clone(f"h_mc_total_stage{i}")
        h3_mc_total.Reset()
        
        for name, h_ptr in mc_hist_ptrs.items():
            h_mc_val = h_ptr.GetValue()
            # Retrieve the pre-cut normalization
            norm_nocut = mc_norms_nocut[name].GetValue() 
            
            if norm_nocut > 0:
                N_fit = fit_results[name]["N"]
                # Scale factor is Fit Yield / Pre-Cut MC Yield
                scale_factor = N_fit / norm_nocut
                h3_mc_total.Add(h_mc_val, scale_factor)
        
        # 3. Calculate Ratio and Inject to C++
        h3_ratio = h3_data.Clone(f"h3_ratio_stage{i}")
        h3_ratio.Divide(h3_mc_total)
        ROOT.gDirectory.Add(h3_ratio)
        
        cpp_code = f"""
        double get_w_stage{i}(double v1, double v2, double v3) {{
            auto h = (TH3D*)gROOT->FindObject("h3_ratio_stage{i}");
            if (!h) return 1.0;
            int b = h->FindBin(v1, v2, v3);
            double w = h->GetBinContent(b);
            
            // Catch 0, empty, nan, or infinity bins
            if (w <= 0.0 || std::isnan(w) || std::isinf(w)) return 1.0;
            
            // Protect against extreme statistical weight spikes 
            if (w > 10.0) return 10.0; 
            
            return w;
        }}
        """
        ROOT.gInterpreter.Declare(cpp_code)
        
        # 4. Multiply current weights by this stage's C++ weight function
        for name in current_weights:
            current_weights[name] = f"({current_weights[name]}) * get_w_stage{i}({vars_3d[0]}, {vars_3d[1]}, {vars_3d[2]})"
            
    # Apply updated weight strings back to the cocktail for final 1D plotting
    for name in cocktail.samples:
        cocktail.samples[name]["w"] = current_weights[name]

# Apply the reweighting
if reweighting_stages:
    apply_sequential_reweighting(data_files, cocktail, fit_results, reweighting_stages)

# --------------------------------------------------
# Plot 
# --------------------------------------------------

output_dir = "25.06.3d_reweighting_plots-full-w"
os.makedirs(output_dir, exist_ok=True)

branch_info = { 
    # B0
    "nTracks": (20, (0.0, 400.0), "nTracks"),
    "b0_PT": (100, (0.0, 40000.0), r"$p_{T}(B^{0})\,[\mathrm{MeV}/c]$"),
    "b0_P": (100, (10000.0, 200000.0), r"$p(B^{0})\,[\mathrm{MeV}/c]$"),
    "b0_ETA": (24, (2.0, 6.0), r"$\eta(B^{0})$"),
    "b0_IP_OWNPV": (50, (0.0, 2.0), r"$\mathrm{IP}(B^{0})\,[\mathrm{mm}]$"),
    "b0_FD_OWNPV": (50, (0.0, 50.0), r"$\mathrm{FD}(B^{0})\,[\mathrm{mm}]$"),
    "b0_FDCHI2_OWNPV": (50, (0.0, 10000.0), r"$\mathrm{FD}\chi^{2}(B^{0})$"),
    "b0_OWNPV_NDOF": (75, (0.0, 300.0), "PV/ndf"),
    "b0_L0Global_TIS": (2, (-0.5, 1.5), "L0Global_TIS"),
    "d0_L0HadronDecision_TOS": (2, (-0.5, 1.5), "L0Hadron_TOS"),

    # Muon
    "mu_PT": (80, (1000.0, 20000.0), r"$p_{T}(\mu)\,[\mathrm{MeV}/c]$"),
    "mu_P": (100, (10000.0, 100000.0), r"$p(\mu)\,[\mathrm{MeV}/c]$"),
    "mu_ETA": (24, (2.0, 6.0), r"$\eta(\mu)$"),

    # D0
    "d0_PT": (100, (2000.0, 15000.0), r"$p_{T}(D^{0})\,[\mathrm{MeV}/c]$"),
    "d0_P": (80, (10000.0, 150000.0), r"$p(D^{0})\,[\mathrm{MeV}/c]$"),
    "d0_ETA": (24, (2.0, 6.0), r"$\eta(D^{0})$"),
    "d0_IP_OWNPV": (50, (0.0, 2.0), r"$\mathrm{IP}(D^{0})\,[\mathrm{mm}]$"),
    "d0_FDCHI2_OWNPV": (50, (0.0, 100000.0), r"$\mathrm{FD}\chi^{2}(D^{0})$"),

    # D*
    "dst_FDCHI2_OWNPV": (50, (0.0, 1000.0), r"$\mathrm{FD}\chi^{2}(D^{*})$"),

    # Kaon
    "k_PT": (80, (900.0, 10000.0), r"$p_{T}(K)\,[\mathrm{MeV}/c]$"),
    "k_P": (100, (0.0, 100000.0), r"$p(K)\,[\mathrm{MeV}/c]$"),
    "k_ETA": (24, (2.0, 6.0), r"$\eta(K)$"),
    "k_IP_OWNPV": (50, (0.0, 5.0), r"$\mathrm{IP}(K)\,[\mathrm{mm}]$"),

    # Pion
    "pi_PT": (80, (900.0, 9000.0), r"$p_{T}(\pi)\,[\mathrm{MeV}/c]$"),
    "pi_P": (100, (0.0, 100000.0), r"$p(\pi)\,[\mathrm{MeV}/c]$"),
    "pi_ETA": (24, (2.0, 6.0), r"$\eta(\pi)$"),
    "pi_IP_OWNPV": (50, (0.0, 5.0), r"$\mathrm{IP}(\pi)\,[\mathrm{mm}]$")

}

# Loop over branches to plot
for branch, (nbins, binrange, xlabel) in branch_info.items():

    fig, (ax, axr) = plt.subplots(
        2, 1,
        figsize=(10, 9),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True
    )

    # Post-fit MC (using the newly sequentially updated weights)
    h_mc, herr_mc, b = cocktail.makePostFitHist(branch, nbins, binrange, fit_results)

    # Data
    h_data, herr_data, _ = makeDataHist(data_files, "DecayTree", branch, nbins, binrange)

    bc = 0.5 * (b[1:] + b[:-1])

    # --- Main plot ---
    #ax.stairs(h_mc, b, label="Post-fit cocktail")
    h_mc_norm = h_mc / np.sum(h_mc) if np.sum(h_mc) > 0 else h_mc
    herr_mc_norm = herr_mc / np.sum(h_mc) if np.sum(h_mc) > 0 else herr_mc
    h_data_norm = h_data / np.sum(h_data) if np.sum(h_data) > 0 else h_data
    herr_data_norm = herr_data / np.sum(h_data) if np.sum(h_data) > 0 else herr_data

    ax.errorbar(bc, h_mc_norm, yerr=herr_mc_norm, fmt="o", markersize=8, capsize=4, label="Post-fit cocktail", color="blue")
  
    ax.errorbar(bc, h_data_norm, yerr=herr_data_norm, fmt="o", color="black", markersize=8, capsize=4, label="Data")

    ax.set_ylabel("Events", fontsize=36)
    ax.set_title(f"$1.\leq w < 1.52$")
    ax.legend(fontsize=34)

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
    axr.errorbar(bc, ratio, yerr=ratio_err, fmt="o", color="black", markersize=8, capsize=4)

    axr.set_xlabel(xlabel, fontsize=56)
    axr.set_ylabel("Data / Model", fontsize=36)
    axr.set_xlim(binrange)
    axr.set_ylim(0.5, 1.5)
    plt.tight_layout()

    output_plot = f"{output_dir}/{branch}_{reweighting_stages[-1]['file_name']}.png"
    plt.savefig(output_plot, dpi=300)
    plt.close(fig)

    print(f"[INFO] Plot saved as {output_plot}")