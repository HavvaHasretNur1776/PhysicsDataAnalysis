#!/usr/bin/env python3

import numpy as np
import pyhf
import matplotlib.pyplot as plt
import json
from redist import modifier
from redist import modifier_hammer
from redist import plot # plotting utilities for hammer components
import sys
import os


### CHOOSE TO FIT REAL DATA OR TOY ####
use_real_data = True

# Merged years data
config_path = f"./fullq2_templates/config_fullq2.json"
data_path = f"./fullq2_templates/RS.dat"
# Load yields 
with open(f'./fullq2_templates/yields/yields.json') as f:
    yield_map = json.load(f)

reader = modifier_hammer.Reader(config_path)
# createFitter produces a fitter object able handle the hammer histograms and produce the templates
fitter = reader.createFitter(False)

# === Template Files CONFIGURATION ===
group_map = {
    "DstMu": [0],
    "DstTau": [1],
#    "DststTau": [2],
    "DstDX": [2],
    "DststHigher": [3],
    "b_dstst": [4], #10413,20413,415,10423,20423,425
    "bs_dstst": [5], #435,10433
    'WSPi':[6], # D* Combinatoric WSMuWSPi NOT removed!
    'WSMu':[7]
}

# === GENERALIZED HELPER ===
def generate_grouped(group_name, params, fitter, mode="toy", ffs_scheme=None):
    """
    Generate combined toy data or template for a group.
    
    :param group_name: name of the group (e.g., 'DstMu')
    :param params: dict of parameters (yields and/or FF deltas)
    :param fitter: fitter object
    :param mode: 'toy' or 'template'
    :param ffs_scheme: optional FF scheme name (e.g., 'BGLVar', 'BLRVar')
    :return: combined template or toy data
    """
    indices = group_map[group_name]
    result = 0

    for idx in indices:
        template = fitter.get_template(idx)

        ffs = {}
        if ffs_scheme == 'BGLVar':
            # Pass FF deltas for BGLVar scheme beside yields
            # results from run1 angular ANA Table 23:
            ffs = {
                "delta_a0": 0.0007, "delta_a1": -0.0007, "delta_a2": 0.0,
                "delta_b0": 0.0, "delta_b1": 0.0017, "delta_b2": 0.0,
                "delta_c1": 0.0004, "delta_c2": -0.0046,"delta_d0": 0.0, "delta_d1": 0.0
            }
        elif ffs_scheme == 'BLRVar':#these we do not use in fit
             ffs = {
                "delta_zt1": 0.0, "delta_ztp": 0.0, "delta_zeta1": 0.0,
                #"delta_chi1": 0.0, "delta_chi2": 0.0, 
                "delta_t1": 0.0, "delta_tp": 0.0, 
                "delta_tau1": 0.0, "delta_tau2": 0.0#, 
                #"delta_eta1": 0.0, "delta_eta2": 0.0, "delta_eta3": 0.0
             }
        # Merge parameters safely — FF deltas only if not already provided
        merged_params = {**ffs, **params}
        if mode == "toy":
#            result += template.generate_toy(**params)
            result += template.generate_toy(**merged_params)
        else:
#            result += template.generate_template(**params)
            result += template.generate_template(**merged_params)
#            print("merged_params:", merged_params)#HERE PARAMS ARE NOT CORRECT
            '''
            Whatever available in yield_map is passed here as: carfeul do not save unused yields to params
            '''
    return result


def generate_toy_dataset(fitter, yield_map):
    """
    Generate combined toy dataset using grouped templates.

    :param fitter: fitter object
    :param yield_map: dictionary with all yields (nDstMu, nDstTau, etc.)
    :return: summed toy dataset
    """
    toy_data = 0

    # Core signal and background components
#    toy_data += generate_grouped("DstMu", {"nDstMu": yield_map["nDstMu"]}, fitter, mode="toy", ffs_scheme="BGLVar")
    toy_data += generate_grouped("DstMu", {"nDstMu": yield_map["nDstMu"]}, fitter, mode="toy")
#    toy_data += generate_grouped("DstTau", {"nDstTau": yield_map["nDstTau"]}, fitter, mode="toy", ffs_scheme="BGLVar")
    toy_data += generate_grouped("DstTau", {"nDstTau": yield_map["nDstTau"]}, fitter, mode="toy")
    #toy_data += generate_grouped("DststTau", {"nDststTau": yield_map["nDststTau"]}, fitter, mode="toy")
    toy_data += generate_grouped("DstDX", {"nDstDX": yield_map["nDstDX"]}, fitter, mode="toy")
    toy_data += generate_grouped("DststHigher", {"nDststHigher": yield_map["nDststHigher"]}, fitter, mode="toy")
    # Grouped D** components (B→D**μν, Bs→D**μν)
    toy_data += generate_grouped("b_dstst", {"nb_dstst": yield_map["nb_dstst"]}, fitter, mode="toy")
    toy_data += generate_grouped("bs_dstst", {"nbs_dstst": yield_map["nbs_dstst"]}, fitter, mode="toy")
    # Wrong-sign pion combinatorial
#    toy_data += generate_grouped("WSPi", {"nWSPi-WSMuWSPi": yield_map["nWSPi-WSMuWSPi"]}, fitter, mode="toy")
    toy_data += generate_grouped("WSPi", {"nWSPi": yield_map["nWSPi"]}, fitter, mode="toy")
    # Wrong-sign muon combinatorial
    toy_data += generate_grouped("WSMu", {"nWSMu": yield_map["nWSMu"]}, fitter, mode="toy")

    return toy_data

# Fitter internally stores toy data whether to fit or not
toy_data = generate_toy_dataset(fitter,yield_map)
fitter.upload_data(toy_data)

# Templates for Null model
def generate_templates(fitter, params_map):#, ffs_scheme=None):
    """
    Generate grouped templates for SM (or NP) under a given FF scheme.

    :param fitter: fitter object
    :param params_map: dict of parameters (e.g. yields)
    :param ffs_scheme: optional string specifying FF scheme, e.g. "BGLVar"
    :return: dict of templates per group
    """
    templates = {}

    # Core components
#    templates["DstMu"] = generate_grouped("DstMu", params_map, fitter, mode="template", ffs_scheme='BGLVar')
    templates["DstMu"] = generate_grouped("DstMu", params_map, fitter, mode="template", ffs_scheme=None)
#    templates["DstTau"] = generate_grouped("DstTau", params_map, fitter, mode="template", ffs_scheme='BGLVar')
    templates["DstTau"] = generate_grouped("DstTau", params_map, fitter, mode="template", ffs_scheme=None)
    #templates["DststTau"] = generate_grouped("DststTau", params_map, fitter, mode="template", ffs_scheme=None)
    templates["DstDX"] = generate_grouped("DstDX", params_map, fitter, mode="template", ffs_scheme=None)
    templates["DststHigher"] = generate_grouped("DststHigher", params_map, fitter, mode="template", ffs_scheme=None)
    # Grouped D** components
    templates["b_dstst"] = generate_grouped("b_dstst", params_map, fitter, mode="template", ffs_scheme=None)
    templates["bs_dstst"] = generate_grouped("bs_dstst", params_map, fitter, mode="template", ffs_scheme=None)
    templates["WSPi"] = generate_grouped("WSPi", params_map, fitter, mode="template", ffs_scheme=None)
    templates["WSMu"] = generate_grouped("WSMu", params_map, fitter, mode="template", ffs_scheme=None)
    ## Debugging
    #print("Params map: ", params_map)
    '''
    Whatever available in yield_map is passed here as: carfeul do not save unused yields to yield_map
    Params map:  {'nDstMu': 258986.4831237528, 'nDstTau': 52942.52710855596, 'nDstDX': 7370.262, 'nb_dstst': 71954.93599433276, 
    'nbs_dstst': 6199.310751779544, 'nDststHigher': 5897.71, 'nWSPi': 1213.0, 'nWSPi-WSMuWSPi': 868.0, 'nWSMu': 85.0}
    '''
    return templates

# Generate NULL theory templates
templates_SM = generate_templates(fitter, yield_map)#, ffs_scheme='BGLVar')

# Define how templates vary in NULL theory (Null dist)
## Wrap them into callables to be used at cmod
def SM_dist(name): return lambda: templates_SM[name]

# Define the alternative model (NP) (the null model is stored in the fitter object)

def NP_dist_DstMu(**kwargs):
#    return generate_grouped("DstMu", kwargs, fitter, mode="template", ffs_scheme="BGLVar")
    return generate_grouped("DstMu", kwargs, fitter, mode="template", ffs_scheme=None)

'''def NP_dist_DstTau(**kwargs):
    return generate_grouped("DstTau", kwargs, fitter, mode="template", ffs_scheme="BGLVar")
'''
def NP_dist_DstTau(**kwargs):
#    print(kwargs.keys())
    f_tau = kwargs["f_Bs_f_tau[1]"]#0.5#kwargs["f_Bs_f_tau[1]"]
    nDstMu = kwargs.get('nDstMu')
    nDstTau = nDstMu * f_tau
#    return generate_grouped("DstTau", {"nDstTau": nDstTau}, fitter, mode="template", ffs_scheme="BGLVar")
    return generate_grouped("DstTau", {"nDstTau": nDstTau}, fitter, mode="template", ffs_scheme=None)

'''def NP_dist_DststTau(**kwargs):
    return generate_grouped("DststTau", kwargs, fitter, mode="template", ffs_scheme=None)
'''
def NP_dist_DstDX(**kwargs):
    return generate_grouped("DstDX", kwargs, fitter, mode="template", ffs_scheme=None)

def NP_dist_DststHigher(**kwargs):
    return generate_grouped("DststHigher", kwargs, fitter, mode="template", ffs_scheme=None)

def NP_dist_b_dstst(**kwargs):
    return generate_grouped("b_dstst", kwargs, fitter, mode="template", ffs_scheme=None)

'''def NP_dist_bs_dstst(**kwargs):
    return generate_grouped("bs_dstst", kwargs, fitter, mode="template", ffs_scheme=None)
'''
def NP_dist_bs_dstst(**kwargs):
    f_Bs = kwargs["f_Bs_f_tau[0]"]
    nb_dstst = kwargs.get('nb_dstst')
    nbs_dstst = nb_dstst * f_Bs 
    return generate_grouped("bs_dstst", {"nbs_dstst": nbs_dstst}, fitter, mode="template", ffs_scheme=None)

def NP_dist_WSPi(**kwargs):
    return generate_grouped("WSPi", kwargs, fitter, mode="template", ffs_scheme=None)

def NP_dist_WSMu(**kwargs):
    return generate_grouped("WSMu", kwargs, fitter, mode="template", ffs_scheme=None)

# Definition of the custom modifiers as the Wilson coefficients and the eigenvectors of the FFs
scale = 0.2
new_params = {
    'nDstMu': {'inits': (yield_map['nDstMu']*scale,), 'bounds': ((0., 1e8),), 'paramset_type': 'unconstrained'},
#    'nDstTau': {'inits': (yield_map['nDstTau']*scale,), 'bounds': ((0., 6e4),), 'paramset_type': 'unconstrained'},
    'nDstDX': {'inits': (yield_map['nDstDX']*scale*2.,), 'bounds': ((0., 5e6),), 'paramset_type': 'unconstrained'},
    'nDststHigher': {'inits': (yield_map['nDststHigher']*3,), 'bounds': ((0., 5e6),), 'paramset_type': 'unconstrained'},
    "nb_dstst": {'inits': (yield_map["nb_dstst"]*scale*0.5,), 'bounds': ((0., 5e6),), 'paramset_type': 'unconstrained'},
#    "nbs_dstst": {'inits': (yield_map["nbs_dstst"]*scale,), 'bounds': ((0., 5e5),), 'paramset_type': 'unconstrained'},    
    'nWSPi': {'inits': (yield_map['nWSPi']*0.2,), 'bounds': ((0., 5e5),), 'paramset_type': 'unconstrained'},
    'nWSMu': {'inits': (yield_map['nWSMu']*2.,), 'bounds': ((0., 5e4),), 'paramset_type': 'unconstrained'},
    "f_Bs_f_tau": {
#            "inits": (0.040,0.04),
            "inits": (0.016,0.0131),
            #N(BsToDsstst)/N(BToDstst): 0.016+/-0.004
            #N(BdToDstTau_TauToMuNuNu)/N(BdToDstMu) 0.0131+/-0.0031
            "bounds": (),
            "cov": [
                [0.004**2, 0.0],
                [0.0, 0.0031**2]
                #[0.010*0.010, 0.0],
                #[0.0, 0.0023*0.0023]#
            ],
            "paramset_type": "constrained_by_normal"
    }
} 
    
# cmod_ created the global object w/ flag "mod${i}" to be used in custom_mod
cmod_sig = modifier_hammer.Modifier_Hammer(new_params, NP_dist_DstMu, SM_dist('DstMu'), name='modS')
cmod_DstTau = modifier_hammer.Modifier_Hammer(new_params, NP_dist_DstTau, SM_dist('DstTau'), name='modT')
cmod_DstDX = modifier_hammer.Modifier_Hammer(new_params, NP_dist_DstDX, SM_dist('DstDX'), name='mod_DstDX')
cmod_DststHigher = modifier_hammer.Modifier_Hammer(new_params, NP_dist_DststHigher, SM_dist('DststHigher'), name='mod_DststHigher')
cmod_b_dstst = modifier_hammer.Modifier_Hammer(new_params, NP_dist_b_dstst, SM_dist('b_dstst'), name='mod_b_dstst')
cmod_bs_dstst = modifier_hammer.Modifier_Hammer(new_params, NP_dist_bs_dstst, SM_dist('bs_dstst'), name='mod_bs_dstst')
cmod_WSPi = modifier_hammer.Modifier_Hammer(new_params, NP_dist_WSPi, SM_dist('WSPi'), name='mod_WSPi')
cmod_WSMu = modifier_hammer.Modifier_Hammer(new_params, NP_dist_WSMu, SM_dist('WSMu'), name='mod_WSMu')

# Map template output to sample names in the spec
sample_data = {
    "nDstMu": templates_SM["DstMu"],
    "nDstTau": templates_SM["DstTau"],
    "nDstDX": templates_SM["DstDX"],
    "nDststHigher": templates_SM["DststHigher"],
    "nb_dstst": templates_SM["b_dstst"],
    "nbs_dstst": templates_SM["bs_dstst"],
    "nWSPi": templates_SM["WSPi"],
    "nWSMu": templates_SM["WSMu"]
}

# Generate samples
samples = [
    {
        "name": name,
        "data": data.tolist(),
        "modifiers": [{"name": "mu", "type": "normfactor", "data": None}]
    }
    for name, data in sample_data.items()
]

# Build spec
spec = {
    "channels": [{"name": "singlechannel","samples": samples}]
}

# Build the model
model = pyhf.Model(spec)

# cutom_mod_ is a dictionary that contains the information about the cmod_ with flag "modS"
## YOU NEED TO ADD THE MODIFIER FOR EACH SAMPLE HERE
custom_mod_sig = {
    "name": "modS_theory",
    "type": "modS",
    "data":{"expr": "modS_weight_fn"}
    }
custom_mod_DstTau = {
    "name": "modT_theory",
    "type": "modT",
    "data":{"expr": "modT_weight_fn"}
    }
'''custom_mod_DststTau = {
    "name": "mod_DststTau_theory",
    "type": "mod_DststTau",
    "data":{"expr": "mod_DststTau_weight_fn"}
    }'''
custom_mod_DstDX = {
    "name": "mod_DstDX_theory",
    "type": "mod_DstDX",
    "data":{"expr": "mod_DstDX_weight_fn"}
    }
custom_mod_DststHigher = {
    "name": "mod_DststHigher_theory",
    "type": "mod_DststHigher",
    "data":{"expr": "mod_DststHigher_weight_fn"}
    }
custom_mod_b_dstst = {
    "name": "mod_b_dstst_theory",
    "type": "mod_b_dstst",
    "data":{"expr": "mod_b_dstst_weight_fn"}
    }
custom_mod_bs_dstst = {
    "name": "mod_bs_dstst_theory",
    "type": "mod_bs_dstst",
    "data":{"expr": "mod_bs_dstst_weight_fn"}
    }
custom_mod_WSPi = {
    "name": "mod_WSPi_theory",
    "type": "mod_WSPi",
    "data":{"expr": "mod_WSPi_weight_fn"}
}
custom_mod_WSMu = {
    "name": "mod_WSMu_theory",
    "type": "mod_WSMu",
    "data":{"expr": "mod_WSMu_weight_fn"}
}
# expaned_pyhf is a dictionary of custom modifiers to be mapped to the pyhf to create pyhf.Model()
expanded_pyhf = {
    **cmod_sig.expanded_pyhf,
    **cmod_DstTau.expanded_pyhf,
    #**cmod_DststTau.expanded_pyhf,
    **cmod_DstDX.expanded_pyhf,
    **cmod_DststHigher.expanded_pyhf,
    **cmod_b_dstst.expanded_pyhf,
    **cmod_bs_dstst.expanded_pyhf,
    **cmod_WSPi.expanded_pyhf,
    **cmod_WSMu.expanded_pyhf
}

# Add components to the Model
model = modifier.add_to_model(model, ['singlechannel'], ['nDstMu'], expanded_pyhf, custom_mod_sig)  
model = modifier.add_to_model(model, ['singlechannel'], ['nDstTau'], expanded_pyhf, custom_mod_DstTau)
model = modifier.add_to_model(model, ['singlechannel'], ['nDstDX'], expanded_pyhf, custom_mod_DstDX) 
model = modifier.add_to_model(model, ['singlechannel'], ['nDststHigher'], expanded_pyhf, custom_mod_DststHigher)
model = modifier.add_to_model(model, ['singlechannel'], ['nb_dstst'], expanded_pyhf, custom_mod_b_dstst)
model = modifier.add_to_model(model, ['singlechannel'], ['nbs_dstst'], expanded_pyhf, custom_mod_bs_dstst)
model = modifier.add_to_model(model, ['singlechannel'], ['nWSPi'], expanded_pyhf, custom_mod_WSPi)
model = modifier.add_to_model(model, ['singlechannel'], ['nWSMu'], expanded_pyhf, custom_mod_WSMu)


# defining the starting values of the d.o.f.
init = model.config.suggested_init()
#init[0] = 0.0 #nDstMu
'''init[1] = 0.0 #DstDX
init[2] = 0.0 #DststHigher
init[3] = 0.0 #b_dstst
'''
#init[4] = 0.0#D* comb
#init[5]= 0.0#WSMu comb
init[8] = 1. #Mu

fixed = model.config.suggested_fixed()
#fixed[0] = True #nDstMu
#fixed[1] = True #nDstDX
#fixed[2] = True #nDststHigher
#fixed[3]=True#nb_dstst
#fixed[4]=True#WSPi D* comb
#fixed[5]=True#WSMu comb
#fixed[6]=True#f_bs
#fixed[7]=True#f_tau
fixed[8] = True #Mu- This has to be here since the normalisation parameter is necessary in pyhf framework even though we do not use it in the fit

# DATA TO BE FITTED
data_obs = np.loadtxt(data_path) if use_real_data else fitter._data.tolist()
#IF YOU HAVE EXTERNAL CONSTRAINTS i.e. gaussian const, you cannot fit ONLY data_obs, you have to add also the auxdata stored in the model.config.auxdata
data = np.concatenate([data_obs, model.config.auxdata]) if use_real_data else fitter._data.tolist() + model.config.auxdata

## FIT 
# pyhf fit and minuit fit
pyhf.set_backend("numpy", pyhf.optimize.scipy_optimizer())
best_fit = pyhf.infer.mle.fit(data, model, fixed_params=fixed, init_pars=init)
pyhf.set_backend("numpy", pyhf.optimize.minuit_optimizer(verbose=True, strategy=2))
best_fit, twice_nll, OptimizeResult = pyhf.infer.mle.fit(
    data,
    model,
    init_pars=best_fit.tolist(),
    fixed_params=fixed,
    return_uncertainties=True,
    return_fitted_val=True,
    return_result_obj=True
    )

# Map parameter names to best-fit values and sigmas
best_fit_dict = modifier.par_dict(model, best_fit)

# separate values and uncertainties from the fit
best_fit_vals = np.array([x[0] for x in best_fit])
best_fit_uncs = np.array([x[1] for x in best_fit])
best_fit_val_dict = modifier.par_dict(model, best_fit_vals)
best_fit_unc_dict = modifier.par_dict(model, best_fit_uncs)

param_names = model.main_model.modifiers_appliers['modS'].param_viewer.allpar_viewer.names
print(f"Fit results in full q2:")
for name in param_names:
    name_str = str(name)

    if "_decorrelated" in name_str:
        continue

    val = best_fit_val_dict.get(name)
    sigma = best_fit_unc_dict.get(name)

    if val is not None and sigma is not None:
        print(f"{name_str} = {val:.6f} ± {sigma:.6f}")

# Getting Gaussian constrained parameters' fitted values and uncertainties in the PHYSICAL space
## When used constrained parameters, the best fit values are in the rotated space
### So we need to rotate them back to the physical space
best_fit_rot_vals = cmod_sig.rotate_pars(best_fit_val_dict)

#### Rotate only values since uncertainties not calculated properly via rotate_pars() method
#best_fit_rot_uncs = cmod_sig.rotate_pars(best_fit_unc_dict)
##### RIGHT WAY TO GET THE UNCERTAINTY ON THE GAUSSIAN CONSTRAINED PARAMETERS IS TO MULTIPLY THE UNROTATED(_decorrelated) UNCERTAINTY BY THE SIGMA OF THE GAUSSIAN CONSTRAINT YOU PUT IN

# Calculate fitted  bs_dstst yields from fitted b_dstst and f_Bs
nb_dstst = best_fit_val_dict['nb_dstst']#no need to use rotated values except for fractions!
fitted_f_Bs = best_fit_rot_vals['f_Bs_f_tau[0]']
fitted_bs_dstst = nb_dstst * fitted_f_Bs 

nDstMu = best_fit_val_dict['nDstMu']
fitted_f_tau = best_fit_rot_vals['f_Bs_f_tau[1]']
fitted_nDstTau = nDstMu * fitted_f_tau

## Error propagation
error_bdstst = best_fit_unc_dict['nb_dstst']
error_f_Bs = best_fit_unc_dict['f_Bs_f_tau_decorrelated[0]']*0.004#use unrotated uncertainty scaled by gaussian constraint sigma that you put
bs_dstst_error = np.sqrt((nb_dstst*error_f_Bs)**2 + (error_bdstst*fitted_f_Bs)**2)

error_nDstMu = best_fit_unc_dict['nDstMu']
error_f_tau = best_fit_unc_dict['f_Bs_f_tau_decorrelated[1]']*0.0031#use unrotated uncertainty scaled by gaussian constraint sigma that you put
fitted_nDstTau_error = np.sqrt( (nDstMu*error_f_tau)**2 + (fitted_f_tau*error_nDstMu)**2 )

# Print derived fitted yields
print(f"nbs_dstst =  {fitted_bs_dstst} +- {bs_dstst_error}")
print(f"nDstTau =  {fitted_nDstTau} +- {fitted_nDstTau_error}")

print(f"f_bs(rotated) =  {fitted_f_Bs} +- {error_f_Bs}")
print(f"f_tau(rotated) =  {fitted_f_tau} +- {error_f_tau}")

# Filter only free parameters (σ > 0)
free_idx = [i for i, s in enumerate(best_fit_uncs) if s > 0]
free_names = [param_names[i] for i in free_idx]

# Get covariance of free parameters only
## MARCO : HESSE probably not handling fixed parameters so covariance for all parameters is always NON-POSITIVE DEFINITE
minuit = OptimizeResult.minuit
cov = minuit.covariance
cov_free = cov[np.ix_(free_idx, free_idx)]

# Eigenvalues for free parameters
eigvals = np.linalg.eigvalsh(cov_free)
if np.all(eigvals > 0):
    print("Covariance matrix of free parameters is positive definite.")
else:
    print("Covariance matrix of free parameters is NOT positive definite.")
    #print("Eigenvalues:", eigvals)

# Correlation matrix for free parameters
corr_free = cov_free / np.outer(np.sqrt(np.diag(cov_free)), np.sqrt(np.diag(cov_free)))
print("Correlation matrix (free parameters only):")
# Pretty print with parameter names
for i, name in enumerate(free_names):
    row = ["{:>8s}".format(name)]
    row += ["{:8.3f}".format(corr_free[i, j]) for j in range(len(free_names))]
    print(" ".join(row))

# Save fitted yields to a JSON file to be used in post-fit
## It can be skipped if not needed
w_edges = [1.000, 1.050, 1.100, 1.150, 1.200, 
           1.250, 1.300, 1.350, 1.400, 1.450, 1.52]

out = {
    "bin": 0,
    "w": {
        "min": w_edges[0],
        "max": w_edges[9],
    },
    "yields": {
        "nDstMu": {
            "val": best_fit_val_dict["nDstMu"],
            "err": best_fit_unc_dict["nDstMu"],
        },
        "nDstTau": {
            "val": fitted_nDstTau,
            "err": fitted_nDstTau_error,
        },
        "nDstDX": {
            "val": best_fit_val_dict["nDstDX"],
            "err": best_fit_unc_dict["nDstDX"],
        },
        "nDststHigher": {
            "val": best_fit_val_dict["nDststHigher"],
            "err": best_fit_unc_dict["nDststHigher"],
        },
        "nb_dstst": {
            "val": best_fit_val_dict["nb_dstst"],
            "err": best_fit_unc_dict["nb_dstst"],
        },
        "nbs_dstst": {
            "val": fitted_bs_dstst,
            "err": bs_dstst_error,
        },
        "nWSPi": {
            "val": best_fit_val_dict["nWSPi"],
            "err": best_fit_unc_dict["nWSPi"],
        },
        "nWSMu": {
            "val": best_fit_val_dict["nWSMu"],
            "err": best_fit_unc_dict["nWSMu"],
        },
    },
}

os.makedirs("fullw_fitted_yields_v2", exist_ok=True)
with open(
    f"fullw_fitted_yields_v2/fittedYields_fullq2.json",
    "w",
) as f:
    json.dump(out, f, indent=2)


#### PLOTTING THE FITTED COMPONENTS AND MODEL ####
best_fit_values = np.array(list(best_fit_rot_vals.values()))# put rotated values
model_counts = model.expected_actualdata(best_fit_values) # expected data counts according to the model
# Define manually parameters names to pass 'best_fit_values' back to fitter to plot components as below
## When you have new parameter or CONSTRAINT you have to update this list accordingly!
param_names_passed_to_fitter=['nDstMu','nDstDX','nDststHigher',"nb_dstst",'nWSPi','nWSMu',"f_Bs_f_tau[0]","f_Bs_f_tau[1]"]#,"f_Bs_f_tau[1]", "f_Bs_f_tau[2]"]#'nWSPi','nWSMu','nWSPi-WSMuWSPi', "f_Bs_f_tau[3]"
bestfit_params = dict(zip(param_names_passed_to_fitter, best_fit_values))

## Generate the fitted components from the best fit values 
components = []
components.append(NP_dist_DstMu(**bestfit_params))
components.append(NP_dist_DstTau(**bestfit_params))
components.append(NP_dist_DstDX(**bestfit_params))
components.append(NP_dist_DststHigher(**bestfit_params))
components.append(NP_dist_b_dstst(**bestfit_params))
components.append(NP_dist_bs_dstst(**bestfit_params))
components.append(NP_dist_WSPi(**bestfit_params))
components.append(NP_dist_WSMu(**bestfit_params))

import mplhep as hep  # LHCb style
# Apply LHCb style
hep.style.use("LHCb2")

# Define components and bin edges
component_names = ["nDstMu", "nDstTau","nDstDX", "nDststHigher",'nb_dstst','nbs_dstst', 'nWSPi','nWSMu']#nWSPi-WSMuWSPi',"nDststTau",'nWSMu','nMisID']
bin_edges = np.linspace(-5.0, 6.0, 21)

bin_centers = (np.array(bin_edges[:-1]) + np.array(bin_edges[1:])) / 2

# Define LaTeX-style labels
component_labels = {
    "nDstMu": r"$B^0 \to D^{*+} \mu^- \nu_\mu$",
    "nDstTau": r"$B^0 \to D^{*+} \tau^- \nu_\tau$",
    "nDststTau": r"$B \to D^{**} \tau^- \nu_\tau$",
    "nDstDX": r"$B \to D^{*}D(X_c \to \mu \nu X)$",
    "nb_dstst": r"$B \to D^{**} \mu^- \nu_\mu$",
    "nbs_dstst": r"$B_s \to D_{s}^{**} \mu^- \nu_\mu$",
    "nDststHigher": r"$B \to D^{**}(Higher) \mu^- \nu_\mu$",
    "nWSMu": "Comb",#"Combinatorial ",
    "nWSPi": r"Comb $D^{*}$",#"Combinatorial D*",
}

# Create output directory for plots
output_dir = "full_q2_fit_projection"
os.makedirs(output_dir, exist_ok=True)
# Define colors
colors = plt.cm.tab10.colors
# Create figure for both subplots
fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]})

# --- Upper subplot: Data, model, and components ---
ax1 = axs[0]
handles = []

# Plot fitted components
for i, (comp_name, comp_counts) in enumerate(zip(component_names, components)):
    line, = ax1.step(bin_edges[:-1], comp_counts, where="post",
                     label=component_labels.get(comp_name, comp_name),
                     linewidth=2, alpha=0.7, color=colors[i % len(colors)])
    handles.append(line)

# Plot total model
ax1.step(bin_edges[:-1], model_counts, where="post", label="Fit", color="red", linestyle="-.", linewidth=2)

# Calculate pulls
pulls = (data_obs - model_counts) / np.sqrt(data_obs) # template uncertainty not propogated nor model uncertainty
ax1.errorbar(bin_centers, data_obs, yerr=np.sqrt(data_obs), fmt="o", color="black", label="Data", markersize=6)

# Styling
ax1.set_ylabel("Counts")
#ax1.set_title("BIN 1")
ax1.set_title(f"full q2 fit")
ax1.legend(fontsize=10)

# --- Lower subplot: Pulls ---
ax2 = axs[1]
ax2.axhline(0, color="black", linestyle="--", linewidth=0.8)
ax2.errorbar(bin_centers, pulls, yerr=np.ones_like(pulls), fmt="o", color="blue")
ax2.set_ylabel("Pulls")
ax2.set_xlabel(r"$m_{miss}^2$ (GeV$^2$)")

# Final layout
plt.tight_layout()
output_plot = f"{output_dir}/fit_projection_RealData_wPulls.png" if use_real_data else f"{output_dir}/fit_projection_ToyData_wPulls.png"

plt.savefig(output_plot, dpi=300)
plt.close()
print(f"Plot saved as {output_plot}")
