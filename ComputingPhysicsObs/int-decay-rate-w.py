import numpy as np
import matplotlib.pyplot as plt
from slophep.Predictions.Observables import BdToDstEllNuPrediction
from slophep.Predictions.SamplingFluctuate import SamplingHelper
from slophep.Predictions.FormFactorsBToV import BdToDstFF
from slophep.utils import setPlotParams
import uproot
import json

# load my unfolded data
my_unf_data_path = "/eos/lhcb/user/h/havva/q2Measurement-BtoDstMuNu/unfolding/unfolded_data_w_iteration_4.root"
my_data = uproot.open(my_unf_data_path)
#bin_edges = np.linspace(1.0, 1.5, 11)
w_edges = [1.000, 1.050, 1.100, 1.150, 1.200,
           1.250, 1.300, 1.350, 1.400, 1.450, 1.520]
w_widths = np.diff(w_edges)
w_data = w_edges[:-1] + w_widths / 2

dGdw_data = my_data["h_unfolded_eff_corrected"].values()
dGdw_err = my_data["h_unfolded_eff_corrected"].errors()

data_integral = np.sum(dGdw_data*w_widths)
dGdw_data_norm = dGdw_data / data_integral
dGdw_err_norm = dGdw_err / data_integral

# Load and Parse HFLAV Data ---
with open('BToDstLNu_Average.json', 'r') as f:
    hflav_json = json.load(f)

# The first measurement is the differential spectra
meas = hflav_json['measurements'][0]
avg_dict = meas['average']

# The first 10 bins correspond to 'w' as per the JSON comment
hflav_bins = np.array(avg_dict['bin'][:10])
hflav_vals = np.array(avg_dict['central'][:10])
hflav_errs = np.array(avg_dict['uncertainty'][:10])

# Calculate bin widths and centers
hflavBin_widths = hflav_bins[:, 1] - hflav_bins[:, 0]
hflavBin_centers = np.mean(hflav_bins, axis=1)

# Normalize to unit area (so sum(height * width) = 1)
# The sum of partial BFs (h_vals) is the total normalization factor
hflav_integral = np.sum(hflav_vals * hflavBin_widths)
hflav_dGdw_norm = hflav_vals / hflav_integral
hflav_dGdw_err_norm = hflav_errs / hflav_integral


# Define constants for B and D* masses (adjust to match your library's values if needed)
mB = 5.279  
mDst = 2.010 

def q2_to_w(q2):
    return (mB**2 + mDst**2 - q2) / (2 * mB * mDst)

def w_to_q2(w):
    return mB**2 + mDst**2 - 2 * mB * mDst * w

def setup_sampler(ff, covdata: str, Nfluct: int):
    pred = BdToDstEllNuPrediction("mu", "mu", ff)
    sampler = SamplingHelper(pred)
    sampler.set_params_from_configfile(covdata)
    sampler.fluctuate(Nfluct)
    return pred, sampler

def get_spectrum_float_w(w_range: list[float], obs, attr: str, fluct):
    res = {"val" : [], "lo" : [], "hi" : []}
    # The Jacobian |dq2/dw| is constant: 2 * mB * mDst
    jacobian = 2 * mB * mDst
    
    for iw in w_range:
        iq2 = w_to_q2(iw)
        o_err = fluct.get_error(attr, [iq2])
        o = getattr(obs, attr)(iq2)
        
        # Convert dG/dq2 to dG/dw
        res["val"].append(o * jacobian)
        res["lo"].append(o_err[0] * jacobian)
        res["hi"].append(o_err[1] * jacobian)
    return res

def plot_spectrum_w(w_range: list[float], 
                    res_list: list[dict], 
                    label_list: list[str], 
                    obs_label: str, 
                    savepath: str = None):
    fig, ax = plt.subplots(1, 1)
    for ires, label in zip(res_list, label_list):
        lines = ax.plot(w_range, ires["val"])
        ax.fill_between(w_range, ires["lo"], ires["hi"], color=lines[0].get_color(), alpha=0.15, label=label)
    
    ax.set(xlabel = r"$w$ (recoil parameter)", ylabel = obs_label)
    ax.legend()
    if savepath:
        fig.savefig(savepath, bbox_inches='tight')
    return fig, ax

### Get the theory predictions and their uncertainties
setPlotParams()

ff_preds = [BdToDstFF.HPQCD,
            BdToDstFF.BGL_FLabMILC, 
            BdToDstFF.BGL_JLQCD,
             BdToDstFF.BGL_Hammer]

labels = ["HPQCD arXiv:2304.03137",
        "Fermilab-MILC arXiv:2105.14019", 
        "JLQCD arXiv:2306.05657",
        "LHCb-CONF-2026-001 (BGL)"]

cov_data = ["slophep/data/BToDstFF_HPQCD_COV_arXiv230403137.json",
    "slophep/data/BToDstFF_BGL_FlabMILC_COV_arXiv210514019.json",
    "slophep/data/BToDstFF_BGL_JLQCD_COV_arXiv230605657.json",
    "BToDstFF_BGLHammer_mymeasurment_systcorr.json"]

Nfluct = 2000
w_range = None
errorbands = []

for iffpred, idat in zip(ff_preds, cov_data):
    ipred, isampler = setup_sampler(iffpred, idat, Nfluct)
    if w_range is None:
        # q2_max corresponds to w=1 (zero recoil)
        # q2_min corresponds to w_max (maximum recoil)
        w_min = q2_to_w(ipred.q2max - 1e-3)
        w_max = q2_to_w(ipred.q2min + 1e-3)
        w_range = np.linspace(w_min, w_max, 500)
    
    errorbands.append(get_spectrum_float_w(w_range, ipred, "dGdq2", isampler))

# Normalize Theory curves so their AREA equals 1
for res in errorbands:
    # Use trapz or simpson to find the area under the theory curve
    theory_integral = np.trapezoid(res["val"], w_range)
    
    res["val"] = np.array(res["val"]) / theory_integral
    res["lo"] = np.array(res["lo"]) / theory_integral
    res["hi"] = np.array(res["hi"]) / theory_integral

fig, ax = plot_spectrum_w(w_range, errorbands, labels, r"$\mathrm{d}\Gamma / \mathrm{d}w$")
# Add my data points with error bars

#ax.errorbar(w_data, dGdw_data_norm, yerr=dGdw_err_norm, fmt='none', color='black', capsize=3, label='Unfolded Data (LHCb)')

## overlay BELLE hflav avg data as well https://gitlab.cern.ch/hflav/hflav-sl-zenodo/-/blob/v1.1.0/BToDstLNu_Average.json
ax.errorbar(hflavBin_centers, hflav_dGdw_norm, 
            yerr=hflav_dGdw_err_norm, 
            xerr=hflavBin_widths/2, 
            fmt='s', color='red', markerfacecolor='none', 
            label='Belle - Belle II HFLAV Avg', capsize=3)
ax.legend()
fig.savefig("bands_dG_dw_belle-avg+run1-NOMYDATA.png", bbox_inches="tight", dpi=300)
