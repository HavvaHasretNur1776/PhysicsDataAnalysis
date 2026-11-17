# from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import json
import scipy
import pyhf
from collections import defaultdict
import fit_setup
from scipy.interpolate import interp1d
from iminuit import Minuit
import time


def zero2epsilon(val, epsilon=0.000001):
    try:
        return np.array([(i if i != 0 else epsilon) for i in val])
    except TypeError:
        return val if val != 0 else epsilon

def get_values_q2_bins(
    bin_idx,
    year,
    extra_label="",
    dim = '1D',
    lumi=5.5,
    jsons_dir=fit_setup.template_path,
    fix_dsttau = False,
    mu_scale=1,
    eff_scale=1,
    toy = False,
    incl = 'all',
    smear = False
):
    num = { k : [] for k in fit_setup.component_templates.keys()}
    num_t = { k : [] for k in fit_setup.component_templates.keys()}
    stat = { k : [] for k in fit_setup.component_templates.keys()}
    uncorr_syst = defaultdict(lambda: np.array([0] * len(fit_setup.mmiss2_bin_edges[f'{dim}_{fit_setup.n_w_bins}'][bin_idx])-1))
    norm_syst = {}
    corr_syst = {"sig": {}, "Dstst": {}, "DstTau": {}, "DstDX": {}, "comb": {}}
    label = f"{year}_q2_{bin_idx}{extra_label}"
    if fix_dsttau:
        label += '_DstTau_fixed'
    out = {}

    with open(fit_setup.template_paths[dim] + f'{year}/' + f'{fit_setup.data_path}{bin_idx}_raw.dat', "r") as f:
        yields = []
        errs = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            v, e = line.split()
            yields.append(float(v))
            errs.append(float(e))
        n_obs = zero2epsilon(yields)
        err_obs = np.array(errs)
    # if toy:
    #     n_obs = np.random.poisson(n_obs)
    #     err_obs = np.sqrt(n_obs)
    scale_obs = np.sum(n_obs)
    #num['obs'] = n_obs
    temp_lab = ''
    w_idx = bin_idx
    if dim == '5D':
        w_idx = bin_idx // 45
    for k,vv in fit_setup.component_templates.items():
            if not smear and k != 'nWSMu':
                temp_lab = '_final'
            else:
                temp_lab = ''
            yields = []
            yields_t = []
            errs = []
            if toy:

                with open(fit_setup.template_paths[dim] + f'{year}/' + f'{vv}{bin_idx}{temp_lab}_raw.dat', "r") as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        v, e = line.split()
                        yields.append(float(v))
                        errs.append(float(e))
                    yields = np.nan_to_num(yields)
                    errs = np.nan_to_num(errs)
                    num[k] = zero2epsilon(yields)
                    stat[k] = np.array(errs)
                    yields_norm = zero2epsilon(yields) / (np.sum(zero2epsilon(yields) ) * np.diff(fit_setup.mmiss2_bin_edges[f'{dim}_{fit_setup.n_w_bins}'][w_idx]))
                    num_t[k] = zero2epsilon(yields_norm)*scale_obs * np.diff(fit_setup.mmiss2_bin_edges[f'{dim}_{fit_setup.n_w_bins}'][w_idx])
            else:
                with open(fit_setup.template_paths[dim] + f'{year}/' + f'{vv}{bin_idx}{temp_lab}_raw.dat', "r") as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        v, e = line.split()
                        yields.append(float(v))
                        errs.append(float(e))
                    num[k] = zero2epsilon(yields)#*scale_obs * np.diff(fit_setup.mmiss2_bin_edges[dim][bin_idx])
                    stat[k] = np.array(errs)#* scale_obs  * np.diff(fit_setup.mmiss2_bin_edges[dim][bin_idx])

    
    if toy:
        print(incl)
        sig_mu = 1
        nn = np.zeros(len(n_obs))
        for k in incl: 
            
            nn += fit_setup.fit_fixed_fracs[k]*num_t[k]  
            sig_mu -= fit_setup.fit_fixed_fracs[k]
            print(k, fit_setup.fit_fixed_fracs[k])
            print(k, num_t[k])
        
        nn += sig_mu*num_t['nDstMu']
        print('nDstMu', sig_mu)
        print('nDstMu', num_t['nDstMu'])
        num['obs'] = np.random.poisson(nn)
    else:
        num['obs'] = n_obs


    # Uncorrelated syst:
    uncorr_syst = { k : np.zeros(len(fit_setup.mmiss2_bin_edges[f'{dim}_{fit_setup.n_w_bins}'][w_idx])-1) for k in fit_setup.component_templates.keys()}
    # Normalisation systs
    norm_syst['nWSMu'] = 0.0
    norm_syst['nWSPi'] = 0.0
    norm_syst["nDstMu"] = 0.0  
    norm_syst["nbDstst"] = 0. 
    norm_syst["nbsDstst"] = 0. 
    norm_syst["nDstst_higher"] = 0. 
    norm_syst["nDstTau"] = 0. 
    norm_syst["nDstDX"] = 0. 


    err = {key: np.sqrt(stat[key] ** 2 + uncorr_syst.get(key, 0) ** 2) for key in stat}

    return num, err, stat, norm_syst, uncorr_syst, corr_syst

def add_channel_q2_bins(
    spec,
    bin_idx,
    year,
    dim = '1D',
    extra_label="",
    syst=True,
    fix_dsttau = False,
    lumi=0.4,
    jsons_dir=fit_setup.template_path,
    mu_scale=1,
    eff_scale=1,
    poi_name=None,
    toy = False,
    incl = 'all',
    smear = False
):

    uncorr_error = "staterror"


    if poi_name is None:
        poi_name = "scale_nDstMu"

    num, err, stat_err, norm_syst, uncorr_syst, corr_syst = get_values_q2_bins(
        bin_idx=bin_idx,
        year = year,
        dim = dim,
        extra_label=extra_label,
        lumi=lumi,
        jsons_dir=jsons_dir,
        fix_dsttau = fix_dsttau,
        mu_scale=mu_scale,
        eff_scale=eff_scale,
        toy = toy,
        incl =incl,
        smear = smear
    )

    channel = {"name": f"w_{bin_idx}", "samples": []}
    samples = {}
    samples["nDstMu"] = {
        "name": f"nDstMu_w_{bin_idx}",
        "data": list(num["nDstMu"]),
        "modifiers": [
            {"name": f'scale_nDstMu_{bin_idx}', "type": "normfactor", "data": None},
         # {"name": 'err_nDstMu', "type": "staterror", "data": err['nDstMu']},
        ],
    }
    if incl == 'all':
        bkg_list = ["nbDstst", "nbsDstst", "nDstst_higher", "nDstTau", "nDstDX", "nWSMu", "nWSPi"]
    else:
        bkg_list = incl
    #bkg_list = ["nbDstst", "nbsDstst", "nDstst_higher", "nDstTau", "nDstDX", "nWSMu", "nWSPi"]
    for bkg in bkg_list:
        if bkg == 'comb':

            samples[bkg] = {
                "name": f"{bkg}_q2_{bin_idx}",
                "data": list(num[bkg]),
                "modifiers": [],
            }
        else:
            if fix_dsttau and bkg == 'nDstTau':
                samples[bkg] = {
                    "name": f"{bkg}_w_{bin_idx}",
                    "data": list(num[bkg]),
                    "modifiers": [],
                }
            else:
                samples[bkg] = {
                    "name": f"{bkg}_w_{bin_idx}",
                    "data": list(num[bkg]),
                    "modifiers": [{"name": f'scale_{bkg}_{bin_idx}', "type": "normfactor", "data": None}]#, {"name": f'err_{bkg}', "type": "staterror", "data": err[bkg]}],
                }

    if syst:
        for name, sample in samples.items():
            # Uncorrelated uncertainties
            sample["modifiers"].append(
                {
                    "name": f"{name}_err_w_{bin_idx}",
                    "type": uncorr_error,
                    "data": list(err[name]),
                }
            )
            # Normalisation uncertainties
            if norm_syst.get(name):
                sample["modifiers"].append(
                    {
                        "name": f"norm_{name}_w_{bin_idx}",
                        "type": "normsys",
                        "data": {"lo": 1 - norm_syst[name], "hi": 1 + norm_syst[name]},
                    }
                )

            for syst_name, values in corr_syst[name].items():
                sample["modifiers"].append(
                    {
                        "name": f"{syst_name}_w_{bin_idx}",
                        "type": "histosys",
                        "data": {
                            "lo_data": list(values[0]),
                            "hi_data": list(values[1]),
                        },
                    }
                )
    else:
        placeholder = 'no stat err for now'

    channel["samples"] = list(samples.values())
    spec["channels"].append(channel)
    
    observation = {"name": f"w_{bin_idx}", "data": list(num["obs"])}
    spec["observations"].append(observation)




def add_channel_q2_bins_simultaneous(
    spec,
    year,
    wbin,
    dim = '1D',
    extra_label="",
    syst=True,
    fix_dsttau = False,
    lumi=0.4,
    jsons_dir=fit_setup.template_path,
    mu_scale=1,
    eff_scale=1,
    poi_name=None,
    toy = False,
    incl = 'all',
    smear = False,
    linear = False,
    quad = False,
):

    uncorr_error = "staterror"
    # uncorr_error = "shapesys"

    if poi_name is None:
        poi_name = "scale_nDstMu"
    if dim == '1D':
        nbins = range(0,fit_setup.n_w_bins)
    else:
        nbins = range(45*wbin, 45*wbin + 45)
        #nbins = fit_setup.nbins_4D['w_Pred'] * fit_setup.nbins_4D['costhetal_reco'] * fit_setup.nbins_4D['costhetad_reco'] * fit_setup.nbins_4D['chi_reco']
    
    for bin_idx in nbins:
        num, err, stat_err, norm_syst, uncorr_syst, corr_syst = get_values_q2_bins(
            bin_idx=bin_idx,
            year = year,
            dim = dim,
            extra_label=extra_label,
            lumi=lumi,
            jsons_dir=jsons_dir,
            fix_dsttau = fix_dsttau,
            mu_scale=mu_scale,
            eff_scale=eff_scale,
            toy = toy,
            incl =incl,
            smear = smear
        )

        channel = {"name": f"w_{bin_idx}", "samples": []}
        samples = {}


        if smear:
            yields_min_sig = []
            yields_max_sig = []
            with open(fit_setup.template_paths[dim] + f'{year}/' + fit_setup.component_templates['nDstMu'] + '{bin_idx}_min_raw.dat', "r") as f:
                for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        v, e = line.split()
                        yields_min_sig.append(float(v))
            with open(fit_setup.template_paths[dim] + f'{year}/' + fit_setup.component_templates['nDstMu']+ f'{bin_idx}_max_raw.dat', "r") as f:
                for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        v, e = line.split()
                        yields_max_sig.append(float(v))
            samples["nDstMu"] = {
                "name": f"nDstMu_w_{bin_idx}",
                "data": list(num["nDstMu"]),
                "modifiers": [
                    {"name": f'scale_nDstMu_{bin_idx}', "type": "normfactor", "data": None},
                                                        {"name": "Resolution",
                                                    "type": "histosys",
                                                    "data": {
                                                        "hi_data": zero2epsilon(yields_max_sig),
                                                        "lo_data": zero2epsilon(yields_min_sig)
                                            }}
                ],
            }
        else:
            samples["nDstMu"] = {
                "name": f"nDstMu_w_{bin_idx}",
                "data": list(num["nDstMu"]),
                "modifiers": [
                    {"name": f'scale_nDstMu_{bin_idx}', "type": "normfactor", "data": None},
                ],
            }
        if incl == 'all':
            bkg_list = ["nbDstst", "nbsDstst", "nDstst_higher", "nDstTau", "nDstDX", "nWSMu", "nWSPi"]
        else:
            bkg_list = incl
        for bkg in bkg_list:
            if bkg == 'nWSMu':
                    samples[bkg] = {
                        "name": f"{bkg}_w_{bin_idx}",
                        "data": list(num[bkg]),
                        "modifiers": [{"name": f'scale_{bkg}_{bin_idx}', "type": "normfactor", "data": None},
                                       #{"name": f'err_{bkg}_{bin_idx}', "type": "staterror", "data": err[bkg]}
                                       
                                       ],
                    }
 
            else:
                if fix_dsttau and bkg == 'nDstTau':
                    samples[bkg] = {
                        "name": f"{bkg}_w_{bin_idx}",
                        "data": list(num[bkg]),
                        "modifiers": [],
                    }
                else:
                    if smear:
                        yields_min = []
                        yields_max = []
                        with open(fit_setup.template_paths[dim] + f'{year}/' + fit_setup.component_templates[bkg] + f'{bin_idx}_min_raw.dat', "r") as f:
                            for line in f:
                                    line = line.strip()
                                    if not line:
                                        continue
                                    v, e = line.split()
                                    yields_min.append(float(v))
                        with open(fit_setup.template_paths[dim] + f'{year}/' + fit_setup.component_templates[bkg] + f'{bin_idx}_max_raw.dat', "r") as f:
                            for line in f:
                                    line = line.strip()
                                    if not line:
                                        continue
                                    v, e = line.split()
                                    yields_max.append(float(v))
                        samples[bkg] = {
                            "name": f"{bkg}_w_{bin_idx}",
                            "data": list(num[bkg]),
                            "modifiers": [{"name": f'scale_{bkg}_{bin_idx}', "type": "normfactor", "data": None},
                                        #    {"name": f'err_{bkg}_{bin_idx}', "type": "staterror", "data": err[bkg]}
                                        {"name": "Resolution",
                                        "type": "histosys",
                                        "data": {
                                            "hi_data": zero2epsilon(yields_max),
                                            "lo_data": zero2epsilon(yields_min)
                                }}
                                        ],
                        }
                    elif quad and bkg == 'nDstDX':
                        yields_min = []
                        yields_max = []
                        with open(fit_setup.template_paths[dim] + f'{year}/' + fit_setup.component_templates[bkg] + f'{bin_idx}_quad_m_raw.dat', "r") as f:
                            for line in f:
                                    line = line.strip()
                                    if not line:
                                        continue
                                    v, e = line.split()
                                    yields_min.append(float(v))
                        with open(fit_setup.template_paths[dim] + f'{year}/' + fit_setup.component_templates[bkg] + f'{bin_idx}_quad_p_raw.dat', "r") as f:
                            for line in f:
                                    line = line.strip()
                                    if not line:
                                        continue
                                    v, e = line.split()
                                    yields_max.append(float(v))
                        samples[bkg] = {
                            "name": f"{bkg}_w_{bin_idx}",
                            "data": list(num[bkg]),
                            "modifiers": [{"name": f'scale_{bkg}_{bin_idx}', "type": "normfactor", "data": None},
                                        #    {"name": f'err_{bkg}_{bin_idx}', "type": "staterror", "data": err[bkg]}
                                        {"name": "Resolution",
                                        "type": "histosys",
                                        "data": {
                                            "hi_data": zero2epsilon(yields_max),
                                            "lo_data": zero2epsilon(yields_min)
                                }}
                                        ],
                        }
                    elif linear and bkg == 'nDstDX':
                        yields_min = []
                        yields_max = []
                        with open(fit_setup.template_paths[dim] + f'{year}/' + fit_setup.component_templates[bkg] + f'{bin_idx}_linear_m_raw.dat', "r") as f:
                            for line in f:
                                    line = line.strip()
                                    if not line:
                                        continue
                                    v, e = line.split()
                                    yields_min.append(float(v))
                        with open(fit_setup.template_paths[dim] + f'{year}/' + fit_setup.component_templates[bkg] + f'{bin_idx}_linear_p_raw.dat', "r") as f:
                            for line in f:
                                    line = line.strip()
                                    if not line:
                                        continue
                                    v, e = line.split()
                                    yields_max.append(float(v))
                        samples[bkg] = {
                            "name": f"{bkg}_w_{bin_idx}",
                            "data": list(num[bkg]),
                            "modifiers": [{"name": f'scale_{bkg}_{bin_idx}', "type": "normfactor", "data": None},
                                        #    {"name": f'err_{bkg}_{bin_idx}', "type": "staterror", "data": err[bkg]}
                                        {"name": "Resolution",
                                        "type": "histosys",
                                        "data": {
                                            "hi_data": zero2epsilon(yields_max),
                                            "lo_data": zero2epsilon(yields_min)
                                }}
                                        ],
                        }
                    else:
                        samples[bkg] = {
                            "name": f"{bkg}_w_{bin_idx}",
                            "data": list(num[bkg]),
                            "modifiers": [{"name": f'scale_{bkg}_{bin_idx}', "type": "normfactor", "data": None},
                                        #    {"name": f'err_{bkg}_{bin_idx}', "type": "staterror", "data": err[bkg]}
                                        ],
                        }

        if syst:
            for name, sample in samples.items():
                # Uncorrelated uncertainties
                sample["modifiers"].append(
                    {
                        "name": f"{name}_err_w_{bin_idx}",
                        "type": uncorr_error,
                        "data": list(err[name]),
                    }
                )
                # Normalisation uncertainties
                if norm_syst.get(name):
                    sample["modifiers"].append(
                        {
                            "name": f"norm_{name}_w_{bin_idx}",
                            "type": "normsys",
                            "data": {"lo": 1 - norm_syst[name], "hi": 1 + norm_syst[name]},
                        }
                    )

                for syst_name, values in corr_syst[name].items():
                    sample["modifiers"].append(
                        {
                            "name": f"{syst_name}_w_{bin_idx}",
                            "type": "histosys",
                            "data": {
                                "lo_data": list(values[0]),
                                "hi_data": list(values[1]),
                            },
                        }
                    )
        else:
            placeholder = 'no stat err for now'

        channel["samples"] = list(samples.values())
        spec["channels"].append(channel)
        
        observation = {"name": f"w_{bin_idx}", "data": list(num["obs"])}
        spec["observations"].append(observation)



def add_measurement_q2_bins(spec, bin_idx,fix_dsttau, dim = '1D', max_range=None, incl = 'all'):
    initial_vals = fit_setup.fit_init_params
    max_range = 1.5 if max_range is None else float(max_range)
    parameters = []

    measurement = {
        "name": 'Measurement',
        "config": {
                "poi": '',
            "parameters": parameters,
        },
    }
    print(measurement)
    spec["measurements"].append(measurement)
    


def add_measurement_q2_bins_simultaneous(spec, nbins, max_range=None, incl = 'all'):
    initial_vals = fit_setup.fit_init_params
    max_range = 1.5 if max_range is None else float(max_range)
    parameters = []
    measurement = {
        "name": 'Measurement',
        "config": {
                "poi": '',
            "parameters": parameters,
        },
    }
    print(measurement)
    spec["measurements"].append(measurement)
    
def get_limit_dict(w, mu_tests=None, cl=0.9, n_points=10, no_systs=False):
    limit_dict = {}
    print('get limits')
    
    bestfit_params, dnll_0 = pyhf.infer.mle.fit(
        w.data(w.model()), w.model(), return_fitted_val=True
    )
    if mu_tests is None:
        mu_tests = auto_mu_tests(w, n_points, list(bestfit_params))


    limit_dict["bestfit_params"] = bestfit_params
    limit_dict["dnll_0"] = dnll_0
    limit_dict["mu_tests"] = mu_tests
    # limit_dict["dnll_tests"] = dnll_tests

    limit_dict["bestfit_params_errors"], limit_dict["bestfit_params_corr"] = get_fit_params_errors(
        w, bestfit_params=list(bestfit_params)
    )
    limit_dict["fitted_yields"] = get_fitted_yields(w, bestfit_params)
    limit_dict["nominal_yields"] = {key : list(val.flatten()) for key, val in zip(w.model().config.samples, w.model().nominal_rates[0])}
    limit_dict["min_nll"] = bestfit_params[w.model().config.poi_index]



    return limit_dict



def auto_mu_tests(w, bin_idx,num, bestfit_params=None):
    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.minuit_optimizer(tolerance=1e-3))
    params = pyhf.infer.mle.fit(
        w.data(w.model()),
        w.model(),
        return_uncertainties=True,
        init_pars=bestfit_params,
    )
    #print(params)
    print(w.model())
    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.scipy_optimizer())
    index = w.model().config.par_names.index(f'scale_nDstMu_{bin_idx}')
    print('index: ',index)
    par = params[index]
    param_dict = {}
    for p,v in zip(w.model().config.par_names,params):
        param_dict[p] = v
    return param_dict

def get_fitted_params(w, bestfit_params=None):
    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.minuit_optimizer(verbose = 1))
    print(pyhf.optimizer.tolerance)
    params, dnll_0 = pyhf.infer.mle.fit(
        w.data(w.model()),
        w.model(),
        return_uncertainties=True,
        return_fitted_val = True,
        init_pars=bestfit_params,
        
    )
    print(dnll_0)
    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.scipy_optimizer())
    index = w.model().config.par_names.index('scale_nDstMu_0')
    print('index: ',index)
    par = params[index]
    param_dict = {}
    for p,v in zip(w.model().config.par_names,params):
        param_dict[p] = v

    limit_dict = {}


    limit_dict["bestfit_params"] = param_dict
    limit_dict["dnll_0"] = dnll_0
  

    limit_dict["fitted_yields"] = get_fitted_yields(w, params[:, 0])
    limit_dict["nominal_yields"] = {key : list(val.flatten()) for key, val in zip(w.model().config.samples, w.model().nominal_rates[0])}

    limit_dict["min_nll"] = params[:, 0]

    return limit_dict


def auto_mu_tests_simultaneous(w, num, dim = '1D', bestfit_params=None):

    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.minuit_optimizer(verbose = True,tolerance=1e-3))
    params, dnll0 = pyhf.infer.mle.fit(
        w.data(w.model()),
        w.model(),
        return_uncertainties=True,
        return_fitted_val = True,
        init_pars=bestfit_params,
    )
    index = w.model().config.par_names.index('scale_nDstMu_0')
    print('index: ',index)
    par = params[index]
    param_dict = {}
    for p,v in zip(w.model().config.par_names,params):
        param_dict[p] = v
    return param_dict

def return_control_yields(year):
    geny = 0
    for s in fit_setup.gen_yields[year]['nDstDX'].keys():
        geny += fit_setup.gen_yields[year]['nDstDX'][s]['up'] + fit_setup.gen_yields[year]['nDstDX'][s]['down']
    fitted_vals = {}
    with open(f'/fit_results_1D/fit_result_control_{year}.txt') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            n, v, e = line.split()
            fitted_vals[n] = [float(v),float(e)]
    print(fitted_vals)
    ddx_yields = []
    ddx_errs = []
    for b in range(5):
        yields_d = []
        errs_d = []
        with open(fit_setup.template_path_control +  f'{year}/'  + f'DstDx_bin{b}_raw.dat', "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                yields_d.append(float(v))
                errs_d.append(float(e))


        ddx_yields.append(np.sum(yields_d))
    eff = np.sum(ddx_yields)/geny
    eff_err = eff * np.sqrt(1/np.sum(ddx_yields) + 1/geny)
    for b in range(5):
        ddx_errs.append(ddx_yields[b] * fitted_vals[f'scale_nDstDX_{b}'][1])
        ddx_yields[b] *= fitted_vals[f'scale_nDstDX_{b}'][0]
    ddx_yields = np.sum(ddx_yields)
    ddx_errs = np.sqrt(np.sum(np.array(ddx_errs)**2))

    return ddx_yields, ddx_errs, eff, eff_err


def constrain_5d_indep(wsp, bin_idx, wbin, toy, year):

    model = wsp.model()
    
    data = wsp.data(model)
    wsmu_constraint = {}
    ddx_constraint = {}


    fit_results_1D = {}
    with open(f'/fit_results_1D/fit_result_1D_{year}.txt') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            n, v, e = line.split()
            if f'_{wbin}' in n:
                fit_results_1D[n] = [float(v),float(e)]
    print(fit_results_1D)
    yield_constraints = {} 

    for k,v in fit_setup.component_templates.items():
        yields = []
        errs = []
        if 'WSMu' in k:
            fn = f'{v}{wbin}_raw.dat'
        else:
            fn = f'{v}{wbin}_final_raw.dat'
        with open(fit_setup.template_paths['1D'] + f'{year}/' + fn , "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                yields.append(float(v))
                errs.append(float(e))
        yield_constraints[k] = [np.sum(yields)* fit_results_1D[f'scale_{k}_{wbin}'][0], np.sum(yields)* fit_results_1D[f'scale_{k}_{wbin}'][1]]
    print(yield_constraints)

    with open(fit_setup.template_paths['1D'] + f'{year}/' + f'data_bin{wbin}_raw.dat', "r") as f:
        yields = []
        errs = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            v, e = line.split()
            yields.append(float(v))
            errs.append(float(e))

    data_wbin = np.sum(yields)

    fracs_ref = {}
    for k,v in fit_setup.component_templates.items():
        fracs_ref[k] = yield_constraints[k][0]/data_wbin
    print("FRACS")
    print(fracs_ref)

    data_yields = {}

    full_templates = {}
    num_bins = range(45*wbin, 45*wbin + 45)
    for binidx in num_bins:
        templates = []
        for k,v in fit_setup.component_templates.items():
            yields = []
            errs = []
            if 'WSMu' in k:
                fn = f'{v}{binidx}_raw.dat'
            else:
                fn = f'{v}{binidx}_final_raw.dat' 
            with open(fit_setup.template_paths['5D'] + f'{year}/' + fn, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    v, e = line.split()
                    yields.append(float(v))
                    errs.append(float(e))
                templates.append(np.array(yields))
                full_templates[f'{k}_{binidx}'] = np.sum(yields)
    sum_yields = {}
    for k in fit_setup.component_templates.keys():
        sum_yields[k] = 0
    for binidx in num_bins:
        for k in fit_setup.component_templates.keys():
            sum_yields[k] += np.nan_to_num(full_templates[f'{k}_{binidx}'])

    yields = []
    errs = []
    with open(fit_setup.template_paths['5D'] + f'{year}/data_bin{bin_idx}_raw.dat', "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            v, e = line.split()
            yields.append(float(v))
            errs.append(float(e))
    data_yields = np.sum(yields)
    constraints = {}
    constraints_toy = {}
    rdst = 0.281
    rdst_e = 0.013
    tautomu = 0.1739
    tautomu_e = 0.0004
    sig_mu = 1
    for k in fit_setup.fit_fixed_fracs.keys():
        sig_mu -= fit_setup.fit_fixed_fracs[k]
    
    for k in fit_setup.component_templates.keys():

        constraints[f'scale_{k}_{bin_idx}'] = [yield_constraints[k][0]/sum_yields[k], yield_constraints[k][1]/ sum_yields[k]]
        if k == 'nDstMu':
            toy_cons = sig_mu*yield_constraints[k][0]/full_templates[f'{k}_{bin_idx}']
        else :
            toy_cons = fit_setup.fit_fixed_fracs[k]*data_yields/full_templates[f'{k}_{bin_idx}']
        if k == 'Dstst_higher' and bin_idx > 179:
            constraints_toy[f'scale_{k}_{bin_idx}'] = [np.nan_to_num(toy_cons), np.nan_to_num(toy_cons)*0.3]
        else:
            constraints_toy[f'scale_{k}_{bin_idx}'] = [np.nan_to_num(toy_cons), np.nan_to_num(toy_cons * constraints[f'scale_{k}_{bin_idx}'][1]/constraints[f'scale_{k}_{bin_idx}'][0])]

    eff_sig = full_templates[f'nDstMu_{bin_idx}']/(fit_setup.gen_yields[year]['nDstMu']['up'] + fit_setup.gen_yields[year]['nDstMu']['down'])
    eff_tau = full_templates[f'nDstTau_{bin_idx}']/(fit_setup.gen_yields[year]['nDstTau']['up'] + fit_setup.gen_yields[year]['nDstTau']['down'])
    eff_sig_e = eff_sig * np.sqrt(1/(full_templates[f'nDstMu_{bin_idx}']) + 1/(fit_setup.gen_yields[year]['nDstMu']['up'] + fit_setup.gen_yields[year]['nDstMu']['down']))
    eff_tau_e = eff_tau * np.sqrt(1/(full_templates[f'nDstTau_{bin_idx}']) + 1/(fit_setup.gen_yields[year]['nDstTau']['up'] + fit_setup.gen_yields[year]['nDstTau']['down']))
    ratio_tau = rdst * tautomu * eff_tau / eff_sig
    ratio_tau_e = ratio_tau * np.sqrt((rdst_e/rdst)**2 + (tautomu_e/tautomu)**2 + (eff_tau_e/eff_tau)**2 + (eff_sig_e/eff_sig)**2)
    constraints[f'MuTau_ratio_{bin_idx}'] = [ratio_tau, ratio_tau_e]
    e_to_r = ratio_tau_e/ratio_tau
    ratio_toy = fit_setup.fit_fixed_fracs['nDstTau']/sig_mu
    constraints_toy[f'MuTau_ratio_{bin_idx}']  = [ratio_toy, e_to_r*ratio_toy]
    print('constraints')
    print(constraints_toy)

    def constrained_loglike(params):

        logL = model.logpdf(params,data)

        
        if toy:         
            p = 0
            p -= 0.5 * ((params[model.config.par_slice(f"scale_nbsDstst_{bin_idx}")][0] -constraints_toy[f'scale_nbsDstst_{bin_idx}'][0]) / (constraints_toy[f'scale_nbsDstst_{bin_idx}'][1])) ** 2
            p -= 0.5 * ((params[model.config.par_slice(f"scale_nbDstst_{bin_idx}")][0] -constraints_toy[f'scale_nbDstst_{bin_idx}'][0]) / (constraints_toy[f'scale_nbDstst_{bin_idx}'][1])) ** 2
            p -= 0.5 * ((params[model.config.par_slice(f"scale_nDstDX_{bin_idx}")][0] -constraints_toy[f'scale_nDstDX_{bin_idx}'][0]) / (constraints_toy[f'scale_nDstDX_{bin_idx}'][1])) ** 2
            p -= 0.5 * ((params[model.config.par_slice(f"scale_nDstst_higher_{bin_idx}")][0] -constraints_toy[f'scale_nDstst_higher_{bin_idx}'][0]) / (constraints_toy[f'scale_nDstst_higher_{bin_idx}'][1])) ** 2
            mu_tau1 =  params[model.config.par_slice(f"scale_nDstTau_{bin_idx}")][0] * full_templates[f"nDstTau_{bin_idx}"]
            mu_tau2 =  params[model.config.par_slice(f"scale_nDstMu_{bin_idx}")][0] * full_templates[f"nDstMu_{bin_idx}"]
        else:
            p = 0

            p -= 0.5 * ((params[model.config.par_slice(f"scale_nbsDstst_{bin_idx}")][0] -constraints[f'scale_nbsDstst_{bin_idx}'][0]) / (constraints[f'scale_nbsDstst_{bin_idx}'][1])) ** 2
            p -= 0.5 * ((params[model.config.par_slice(f"scale_nbDstst_{bin_idx}")][0] -constraints[f'scale_nbDstst_{bin_idx}'][0]) / (constraints[f'scale_nbDstst_{bin_idx}'][1])) ** 2
            p -= 0.5 * ((params[model.config.par_slice(f"scale_nDstDX_{bin_idx}")][0] -constraints[f'scale_nDstDX_{bin_idx}'][0]) / (constraints[f'scale_nDstDX_{bin_idx}'][1])) ** 2
            mu_tau1 =  params[model.config.par_slice(f"scale_nDstTau_{bin_idx}")][0] * full_templates[f"nDstTau_{bin_idx}"]
            mu_tau2 =  params[model.config.par_slice(f"scale_nDstMu_{bin_idx}")][0] * full_templates[f"nDstMu_{bin_idx}"]
            p -= 0.5 * (((mu_tau1/mu_tau2) - constraints[f'MuTau_ratio_{bin_idx}'][0]) / (constraints[f'MuTau_ratio_{bin_idx}'][1]))**2


        return float(np.sum(logL) + p)

    init_pars = model.config.suggested_init()


    def fcn(*pars):
        return -constrained_loglike(pars)
    name = model.config.par_order

    m = Minuit(fcn, *init_pars, name=name)

    for i in name:
        m.limits[i] = [0.0, 10.0]

    m.errordef = 0.5 


    m.migrad()
    #cov = np.array(m.covariance)
    m.hesse()
    param_dict = {}
    params = []

    params_ref = []
    for k in model.config.par_order:
        params_ref.append(m.values[k])
        param_dict[k] = [m.values[k], m.errors[k]]
        params.append(np.array([m.values[k], m.errors[k]]))
    params = np.array(params)
    dnll_0 = m.fval
    print(dnll_0)
    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.scipy_optimizer())


    limit_dict = {}


    limit_dict["bestfit_params"] = param_dict
    limit_dict["dnll_0"] = dnll_0
  

    limit_dict["fitted_yields"] = get_fitted_yields(wsp, params[:, 0])
    limit_dict["nominal_yields"] = {key : list(val.flatten()) for key, val in zip(wsp.model().config.samples, wsp.model().nominal_rates[0])}

    #limit_dict["min_nll"] = bestfit_params[w.model().config.poi_index]
    limit_dict["min_nll"] = params[:, 0]
    print('Param dict: ',param_dict)
    print('params ref: ', params_ref)

    if toy:
        print('DATA:')
        print(np.sum(data))
        print(np.sum(model.spec['channels'][0]['samples'][0]['data']) * param_dict[f'scale_nDstMu_{bin_idx}'][0]/np.sum(data))
        sigfrac = np.sum(model.spec['channels'][0]['samples'][0]['data']) * param_dict[f'scale_nDstMu_{bin_idx}'][0]/np.sum(data)
        sigerr = np.sum(model.spec['channels'][0]['samples'][0]['data']) * param_dict[f'scale_nDstMu_{bin_idx}'][1]/np.sum(data)
        sig_mu = 1
        for k, v in fit_setup.fit_fixed_fracs.items():   
            sig_mu -= v
        print('nDstMu', sig_mu)
        limit_dict = [sig_mu, sigfrac, sigerr]#, cov, m.parameters]
    return limit_dict
  
def constrain_bs_simultaneous(wsp, num_bins, toy, year, smear):
    model = wsp.model()
    
    data = wsp.data(model)
    wsmu_constraint = {}
    ddx_constraint = {}
    
    bs_toty = []
    b_toty = []
    sig_toty = []
    tau_toty = []
    yields_d = []
    yields_w = []
    yields_ddx = []
    
    for bin_number in range(num_bins):
        yw = 0
        yd = 0
        ybs = 0
        ybd = 0
        ysig = 0
        ytau = 0
        yddx = 0
        with open(fit_setup.template_paths['1D'] + f'{year}/' + f'data_WSMu_bin{bin_number}_raw.dat', "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                yw += float(v)
        with open(fit_setup.template_paths['1D'] + f'{year}/' + f'DstDx_bin{bin_number}_raw.dat', "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                yddx += float(v)
        with open(fit_setup.template_paths['1D'] + f'{year}/' + f'b_Dstst_bin{bin_number}_raw.dat', "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                ybd += float(v)
        with open(fit_setup.template_paths['1D'] + f'{year}/' + f'bs_Dstst_bin{bin_number}_raw.dat', "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                ybs += float(v)
        with open(fit_setup.template_paths['1D'] + f'{year}/' + f'b0_DstTau_bin{bin_number}_raw.dat', "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                ytau += float(v)
        with open(fit_setup.template_paths['1D'] + f'{year}/' + f'signal_bin{bin_number}_raw.dat', "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                ysig+= float(v)
        with open(fit_setup.template_paths['1D'] + f'{year}/' + f'{fit_setup.data_path}{bin_number}_raw.dat', "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                yd +=float(v)
        wsmu_constraint[f'ratio_{bin_number}'] = yw/yd
        ddx_constraint[f'ratio_{bin_number}'] = yddx/yd
        yields_w.append(yw)
        yields_d.append(yd)
        yields_ddx.append(yddx)
        bs_toty.append(ybs)
        b_toty.append(ybd)
        sig_toty.append(ysig)
        tau_toty.append(ytau)
    fracs_tau = []
    for i in range(5):
        fracs_tau.append(tau_toty[i]/np.sum(tau_toty))
    print('ddx constraint: ',ddx_constraint)
    print(wsmu_constraint)
    print('len data ', len(data))
    if toy:
        bs_f = np.array([fit_setup.fit_fixed_fracs['nbsDstst'] * k for k in yields_d])/(np.sum(yields_d) * fit_setup.fit_fixed_fracs['nbsDstst'])
    else:
        bs_f = bs_toty/np.sum(bs_toty)
        bs_f_e  = []
        for bsf in bs_f:
            bs_f_e.append(np.sqrt((bsf*(1-bsf))/np.sum(bs_toty)))
    bs_gen = 0
    for ssg, ssgy in fit_setup.gen_yields[year]['nbsDstst'].items():
        bs_gen += ssgy['up'] + ssgy['down']
    b_gen = 0
    for ssg, ssgy in fit_setup.gen_yields[year]['nbDstst'].items():
        b_gen += ssgy['up'] + ssgy['down']
    bs_eff = np.sum(bs_toty)/bs_gen
    bs_eff_e = bs_eff * np.sqrt(1/np.sum(bs_toty) + 1/bs_gen)
    b_eff = np.sum(b_toty)/b_gen
    b_eff_e = b_eff * np.sqrt(1/np.sum(b_toty) + 1/b_gen)
    bs_b_ratio = 0.16 * bs_eff / b_eff
    bs_b_e = bs_b_ratio * np.sqrt((0.04/0.16)**2 + (b_eff_e/b_eff)**2 + (bs_eff_e/bs_eff)**2)


    control_yields, control_err, control_eff, control_eff_err = return_control_yields(year)
    ddx_gen = 0
    for ssg, ssgy in fit_setup.gen_yields[year]['nDstDX'].items():
        ddx_gen += ssgy['up'] + ssgy['down']
    sr_eff = np.sum(yields_ddx)/ddx_gen
    sr_eff_e = sr_eff * np.sqrt(1/np.sum(yields_ddx) + 1/ddx_gen)
    exp_ddx = control_yields * sr_eff / control_eff
    exp_ddx_e = exp_ddx * np.sqrt((control_eff/control_yields)**2 + (control_eff_err/control_eff)**2 + (sr_eff/sr_eff_e)**2)

    rdst = 0.281
    rdst_e = 0.013
    tautomu = 0.1739
    tautomu_e = 0.0004
    eff_sig = np.array(sig_toty)/(fit_setup.gen_yields[year]['nDstMu']['up'] + fit_setup.gen_yields[year]['nDstMu']['down'])
    eff_tau = np.array(tau_toty)/(fit_setup.gen_yields[year]['nDstTau']['up'] + fit_setup.gen_yields[year]['nDstTau']['down'])
    eff_sig_e = eff_sig * np.sqrt(1/(np.array(sig_toty)) + 1/(fit_setup.gen_yields[year]['nDstMu']['up'] + fit_setup.gen_yields[year]['nDstMu']['down']))
    eff_tau_e = eff_tau * np.sqrt(1/(np.array(tau_toty)) + 1/(fit_setup.gen_yields[year]['nDstTau']['up'] + fit_setup.gen_yields[year]['nDstTau']['down']))
    ratio_tau = rdst * tautomu * eff_tau / eff_sig
    ratio_tau_e = ratio_tau * np.sqrt((rdst_e/rdst)**2 + (tautomu_e/tautomu)**2 + (eff_tau_e/eff_tau)**2 + (eff_sig_e/eff_sig)**2)

    
    def constrained_loglike(params):

        logL = model.logpdf(params,data)
        mu1 = 0
        mu2 = 0
        mu_tau1 = 0
        mu_tau2 = 0
        sum_bs = 0
        sum_tau = 0
        sum_ddx = 0
        
        if toy:
            
            p = 0
            for bin_number in range(num_bins):
                r = fit_setup.fit_fixed_fracs['nWSMu']
                p -= 0.5 * (( params[model.config.par_slice(f"scale_nWSMu_{bin_number}")][0]*yields_w[bin_number]/yields_d[bin_number] - r) / (r_e)) ** 2

                p -= 0.5 * (( params[model.config.par_slice(f"scale_nDstDX_{bin_number}")][0]*yields_ddx[bin_number]/yields_d[bin_number] - ref_dstdx) / (ref_dstdx_e)) ** 2
            for bin_number in range(num_bins):
                mu1 += params[model.config.par_slice(f"scale_nbsDstst_{bin_number}")][0] * bs_toty[bin_number]
                sum_bs += params[model.config.par_slice(f"scale_nbsDstst_{bin_number}")][0]  * bs_toty[bin_number]
                sum_tau += params[model.config.par_slice(f"scale_nDstTau_{bin_number}")][0]  * tau_toty[bin_number]
                mu2 += params[model.config.par_slice(f"scale_nbDstst_{bin_number}")][0] * b_toty[bin_number]
                mu_tau1 +=  params[model.config.par_slice(f"scale_nDstTau_{bin_number}")][0] * tau_toty[bin_number]
                mu_tau2 +=  params[model.config.par_slice(f"scale_nDstMu_{bin_number}")][0] * sig_toty[bin_number]
            p -= 0.5 * (((mu1/mu2) - ratio) / ratio_e) ** 2

            for bin_number in range(num_bins):
                p -= 0.5 * (((params[model.config.par_slice(f"scale_nbsDstst_{bin_number}")][0]*bs_toty[bin_number]/sum_bs) - bs_f[bin_number]) / (bs_f[bin_number]*bs_err_fac))**2
                p -= 0.5 * (((params[model.config.par_slice(f"scale_nDstTau_{bin_number}")][0]*tau_toty[bin_number]) - tau_fix[bin_number]) / (tau_fix[bin_number]*tau_fix_err_fac))**2
    
        else:


            p = 0

            for bin_number in range(num_bins):
                mu1 += params[model.config.par_slice(f"scale_nbsDstst_{bin_number}")][0] * bs_toty[bin_number]
                sum_bs += params[model.config.par_slice(f"scale_nbsDstst_{bin_number}")][0]  * bs_toty[bin_number]
                #sum_tau += params[model.config.par_slice(f"scale_nDstTau_{bin_number}")][0]  * tau_toty[bin_number]
                sum_ddx += params[model.config.par_slice(f"scale_nDstDX_{bin_number}")][0]  * yields_ddx[bin_number]
                mu2 += params[model.config.par_slice(f"scale_nbDstst_{bin_number}")][0] * b_toty[bin_number]
                mu_tau1 =  params[model.config.par_slice(f"scale_nDstTau_{bin_number}")][0] * tau_toty[bin_number]
                mu_tau2 =  params[model.config.par_slice(f"scale_nDstMu_{bin_number}")][0] * sig_toty[bin_number]
                p -= 0.5 * (((mu_tau1/mu_tau2) - ratio_tau[bin_number]) / (ratio_tau_e[bin_number]))**2

            p -= 0.5 * (((mu1/mu2) - bs_b_ratio) / (bs_b_e)) ** 2
            p -= 0.5 * ((sum_ddx - exp_ddx) / (exp_ddx_e)) ** 2
 
            for bin_number in range(num_bins):
                p -= 0.5 * (((params[model.config.par_slice(f"scale_nbsDstst_{bin_number}")][0]*bs_toty[bin_number]/sum_bs) - bs_f[bin_number]) / (bs_f_e[bin_number]))**2


        return float(np.sum(logL) + p)
    init_pars = model.config.suggested_init()
    print('par_order: ', model.config.par_order)
    print('init pars: ',init_pars)
    npars = len(init_pars)
    bounds = model.config.suggested_bounds() 

    init_pars = np.array(init_pars)
    #print(bounds)
    # print(len(in_pars))
    def fcn(*pars):
        return -constrained_loglike(pars)
    name = model.config.par_order
    m = Minuit(fcn, *init_pars, name=name)
    for i in name:
        if "WSMu" in i:
            m.limits[i] = [0.0, 1000.0]
        else:
            m.limits[i] = [0.0, 10.0]
    m.errordef = 0.5 
    #m.tol = 0.001
    m.migrad()
    m.hesse()
    param_dict = {}
    params = []
    params_ref = []
    for k in model.config.par_order:
        params_ref.append(m.values[k])
        param_dict[k] = [m.values[k], m.errors[k]]
        params.append(np.array([m.values[k], m.errors[k]]))
    params = np.array(params)
    dnll_0 = m.fval
    print(dnll_0)
    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.scipy_optimizer())


    limit_dict = {}


    limit_dict["bestfit_params"] = param_dict
    limit_dict["dnll_0"] = dnll_0
  

    limit_dict["fitted_yields"] = get_fitted_yields(wsp, params[:, 0])
    limit_dict["nominal_yields"] = {key : list(val.flatten()) for key, val in zip(wsp.model().config.samples, wsp.model().nominal_rates[0])}

    #limit_dict["min_nll"] = bestfit_params[w.model().config.poi_index]
    limit_dict["min_nll"] = params[:, 0]
    print('Param dict: ',param_dict)
    print('params ref: ', params_ref)

    with open(fit_setup.template_paths['1D'] + f'{year}/' + f'data_WSMu_bin{bin_number}_raw.dat', "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            v, e = line.split()
            yw += float(v)
    if smear:
        def set_resolution(model, base_pars, value):
            pars_ = np.array(base_pars, copy=True)
            idx = model.config.par_order.index("Resolution")
            pars_[idx] = value
            return pars_
        fitted   = model.expected_actualdata(set_resolution(model, params_ref, limit_dict['bestfit_params']['Resolution'][0]))
        print(fitted)
        print(len(fitted))
        sample_contributions = {}

        for sample in model.config.samples:
            pars = params_ref.copy()
            for i, parname in enumerate(model.config.par_order):
                #print(parname)
                if parname not in ['scale_' + sample.replace('_w',''), 'Resolution']:
                        pars[i] = 0.0  # remove other samples

            pred = model.expected_actualdata(pars)
            sample_contributions[sample] = pred

        grouped = defaultdict(list)

        for name, arr in sample_contributions.items():
            idx = int(name.split("_")[-1])  # extract 0–4
            grouped[idx].append((name.replace(f'_w_{idx}',''), arr))
        n_per_group = 20
        print('CHANNELS')
        print(wsp['channels'][0])

        def return_template(name, binno):
            channel = wsp['channels'][binno]
            sample = next(s['data'] for s in channel["samples"] if s["name"] == f'{name}_w_{binno}')
            return np.array(sample)

        for idx in sorted(grouped):
            yields = []
            errs = []
            bin_edges = fit_setup.mmiss2_bin_edges[f'1D_{fit_setup.n_w_bins}'][idx]
            bin_centers = (np.array(bin_edges[:-1]) + np.array(bin_edges[1:])) / 2
            bin_widths = np.diff(bin_edges)
            with open(fit_setup.template_paths['1D'] + f'{year}/' + f'{fit_setup.data_path}{idx}_raw.dat', "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    v, e = line.split()
                    yields.append(float(v))
                    errs.append(float(e))
            data= np.array(yields)
            data_err = np.array(errs)
            start = idx * n_per_group
            end = (idx + 1) * n_per_group

            bins = np.arange(n_per_group)
            bottom = np.zeros(n_per_group)
            #full_err = np.zeros(n_per_group)
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True ,gridspec_kw={'height_ratios': [3, 1]})
            print('GROUPED')
            print(grouped[idx])
            vals = {}
            for name, arr in reversed(grouped[idx]):
                if 'DstMu' in name:

                    print(name, arr)
                vals[name] = arr[start:end]
            for name in ["nbDstst", "nbsDstst", "nDstst_higher", "nDstTau", "nDstDX", 'nWSMu', "nDstMu"]:
                ax1.bar(bin_centers, vals[name], width=bin_widths,bottom=bottom, color = fit_setup.component_colors[name],label=fit_setup.component_labels[name]+F': {np.sum(vals[name]):.0f}')
                bottom += vals[name]
            ax1.errorbar(bin_centers, data, yerr = data_err, xerr = bin_widths/2, color = 'black', label = 'data' + f': {np.sum(data):.0f}', fmt = 'o',  ms = 2, capsize = 1)
            ax1.set_ylim(0,1.1*bottom.max())
            ax1.set_xlabel("Bin")
            ax1.set_ylabel("Events")
            ax1.legend(reverse = True)
            
            full_err =np.zeros(20)
            for name in ["nbDstst", "nbsDstst", "nDstst_higher", "nDstTau", "nDstDX", 'nWSMu', "nDstMu"]:
                full_err += (param_dict[f"scale_{name}_{idx}"][1] * return_template(name,idx))**2
            full_err = np.sqrt(full_err)
            for i, err in zip(range(len(bin_centers)), full_err):
                x = bin_edges[i]
                height = bottom[i]
                ax1.add_patch(plt.Rectangle((x, height-np.sqrt(err)/2), bin_widths[i], np.sqrt(err), fill=True, facecolor = 'darkgrey', alpha = 0.4,
                                        edgecolor='black', hatch='////', lw = 0))
            pulls = (bottom- data)/np.sqrt( full_err**2 + data_err**2)
            #pulls = (bot- data)/np.sqrt(data_err**2)

            ax2.axhline(0, color="black", linestyle="--", linewidth=0.8)
            ax2.errorbar(bin_centers, pulls, yerr=np.ones_like(pulls), fmt="o", color="black")
            ax2.set_ylim(-5,5)
            ax2.set_ylabel("Pulls")
            plt.savefig(f'check_smear_bin{idx}_{year}.pdf')
            plt.clf()
            for name, arr in reversed(grouped[idx]):
                if 'WSMu' not in name:
                    yyy = []
                    with open(fit_setup.template_paths['1D'] + f'{year}/' + f'{fit_setup.component_templates[name]}{idx}_final_raw.dat', "r") as f:
                        for line in f:
                            line = line.strip()
                            if not line:
                                continue
                            v, e = line.split()
                            yyy.append(float(v))

                        yyy = np.array(yyy)
                    yyy = yyy/np.sum(yyy)
                    vals = arr[start:end]
                    vals = np.array(vals)/np.sum(vals)
                    fig, ax1 = plt.subplots(1, 1, figsize=(10, 8))

                    ax1.bar(bin_centers, vals, width=bin_widths, color = fit_setup.component_colors[name],label=fit_setup.component_labels[name] + ' smearing fitted')
                    ax1.step(np.append(bin_edges, bin_edges[-1])[:-1], np.append(yyy,yyy[-1]), where="post",
                                label= f'smearing applied',
                                linewidth=1, color='black') 
                    ax1.set_xlabel("Bin")
                    ax1.set_ylabel("A.U")
                    ax1.legend()
                    #plt.yscale('log')
                    plt.savefig(f'smearing_checks/check_smear_{name}_bin{idx}_{year}.pdf')
                    plt.clf()

    return limit_dict


def constrain_5d_simultaneous(wsp, wbin, toy, year):

    model = wsp.model()
    
    data = wsp.data(model)
    wsmu_constraint = {}
    ddx_constraint = {}


    fit_results_1D = {}
    with open(f'/fit_results_1D/fit_result_1D_{year}.txt') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            n, v, e = line.split()
            if f'_{wbin}' in n:
                fit_results_1D[n] = [float(v),float(e)]
    print(fit_results_1D)
    yield_constraints = {} 

    for k,v in fit_setup.component_templates.items():
        yields = []
        errs = []
        if 'WSMu' in k:
            fn = f'{v}{wbin}_raw.dat'
        else:
            fn = f'{v}{wbin}_final_raw.dat'
        with open(fit_setup.template_paths['1D'] + f'{year}/' + fn , "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                yields.append(float(v))
                errs.append(float(e))
        yield_constraints[k] = [np.sum(yields)* fit_results_1D[f'scale_{k}_{wbin}'][0], np.sum(yields)* fit_results_1D[f'scale_{k}_{wbin}'][1]]
    print(yield_constraints)

    with open(fit_setup.template_paths['1D'] + f'{year}/' + f'data_bin{wbin}_raw.dat', "r") as f:
        yields = []
        errs = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            v, e = line.split()
            yields.append(float(v))
            errs.append(float(e))

    data_wbin = np.sum(yields)

    fracs_ref = {}
    for k,v in fit_setup.component_templates.items():
        fracs_ref[k] = yield_constraints[k][0]/data_wbin
    print("FRACS")
    print(fracs_ref)
    data_obs = []
    bs_toty = []
    b_toty = []
    sig_toty = []
    tau_toty = []
    yields_d = []
    yields_w = []
    yields_ddx = []
    num_bins = range(45*wbin, 45*wbin + 45)

    indep_res = {}
    data_yields = {}

    for bin_idx in num_bins:
        with open(f'results_5D/fit_result_5D_{bin_idx}_{year}.txt') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                n, v, e = line.split()
 
                indep_res[n] = [float(v),float(e)]

        out = {}
        with open(fit_setup.template_paths['5D'] + f'{year}/' + f'data_bin{bin_idx}_raw.dat', "r") as f:
            yields = []
            errs = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                yields.append(float(v))
                errs.append(float(e))
            n_obs = np.array(yields)
            err_obs = np.array(errs)
        data_obs.append(n_obs)
        data_yields[bin_idx] = n_obs

    full_templates = {}
    for bin_idx in num_bins:
        templates = []
        for k,v in fit_setup.component_templates.items():
            yields = []
            errs = []
            if 'WSMu' in k:
                fn = f'{v}{bin_idx}_raw.dat'
            else:
                fn = f'{v}{bin_idx}_final_raw.dat' 
            with open(fit_setup.template_paths['5D'] + f'{year}/' + fn, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    v, e = line.split()
                    yields.append(float(v))
                    errs.append(float(e))
                templates.append(np.array(yields))
                full_templates[f'{k}_{bin_idx}'] = np.sum(yields)
    sum_yields = {}
    for k in fit_setup.component_templates.keys():
        sum_yields[k] = 0
    for bin_idx in num_bins:
        for k in fit_setup.component_templates.keys():
            sum_yields[k] += full_templates[f'{k}_{bin_idx}']
    upper_bounds = {}
    constraints = {}
    rdst = 0.281
    rdst_e = 0.013
    tautomu = 0.1739
    tautomu_e = 0.0004
    for bin_idx in num_bins:
        for k in fit_setup.component_templates.keys():
            upper_bounds[f'scale_{k}_{bin_idx}'] = 1.1 * np.sum(data_yields[bin_idx])/np.sum(full_templates[f'{k}_{bin_idx}'])
            constraints[f'scale_{k}_{bin_idx}'] = [yield_constraints[k][0]/sum_yields[k], yield_constraints[k][1]/ sum_yields[k]]
    for bin_idx in num_bins:
        eff_sig = full_templates[f'nDstMu_{bin_idx}']/(fit_setup.gen_yields[year]['nDstMu']['up'] + fit_setup.gen_yields[year]['nDstMu']['down'])
        eff_tau = full_templates[f'nDstTau_{bin_idx}']/(fit_setup.gen_yields[year]['nDstTau']['up'] + fit_setup.gen_yields[year]['nDstTau']['down'])
        eff_sig_e = eff_sig * np.sqrt(1/(full_templates[f'nDstMu_{bin_idx}']) + 1/(fit_setup.gen_yields[year]['nDstMu']['up'] + fit_setup.gen_yields[year]['nDstMu']['down']))
        eff_tau_e = eff_tau * np.sqrt(1/(full_templates[f'nDstTau_{bin_idx}']) + 1/(fit_setup.gen_yields[year]['nDstTau']['up'] + fit_setup.gen_yields[year]['nDstTau']['down']))
        ratio_tau = rdst * tautomu * eff_tau / eff_sig
        ratio_tau_e = ratio_tau * np.sqrt((rdst_e/rdst)**2 + (tautomu_e/tautomu)**2 + (eff_tau_e/eff_tau)**2 + (eff_sig_e/eff_sig)**2)
        constraints[f'MuTau_ratio_{bin_idx}'] = [ratio_tau, ratio_tau_e]
    print('constraints')
    print(constraints)

    def constrained_loglike(params):

        logL = model.logpdf(params,data)
        sum_bs = 0
        sum_bd = 0
        sum_ddx = 0
        
        if toy:
            
            p = 0
            for bin_number in range(num_bins):
                r = fit_setup.fit_fixed_fracs['nWSMu']
                p -= 0.5 * (( params[model.config.par_slice(f"scale_nWSMu_{bin_number}")][0]*yields_w[bin_number]/yields_d[bin_number] - r) / (r_e)) ** 2

                p -= 0.5 * (( params[model.config.par_slice(f"scale_nDstDX_{bin_number}")][0]*yields_ddx[bin_number]/yields_d[bin_number] - ref_dstdx) / (ref_dstdx_e)) ** 2
            for bin_number in range(num_bins):
                mu1 += params[model.config.par_slice(f"scale_nbsDstst_{bin_number}")][0] * bs_toty[bin_number]
                sum_bs += params[model.config.par_slice(f"scale_nbsDstst_{bin_number}")][0]  * bs_toty[bin_number]
                sum_tau += params[model.config.par_slice(f"scale_nDstTau_{bin_number}")][0]  * tau_toty[bin_number]
                mu2 += params[model.config.par_slice(f"scale_nbDstst_{bin_number}")][0] * b_toty[bin_number]
                mu_tau1 +=  params[model.config.par_slice(f"scale_nDstTau_{bin_number}")][0] * tau_toty[bin_number]
                mu_tau2 +=  params[model.config.par_slice(f"scale_nDstMu_{bin_number}")][0] * sig_toty[bin_number]
            p -= 0.5 * (((mu1/mu2) - ratio) / ratio_e) ** 2
            for bin_number in range(num_bins):
                p -= 0.5 * (((params[model.config.par_slice(f"scale_nbsDstst_{bin_number}")][0]*bs_toty[bin_number]/sum_bs) - bs_f[bin_number]) / (bs_f[bin_number]*bs_err_fac))**2
                p -= 0.5 * (((params[model.config.par_slice(f"scale_nDstTau_{bin_number}")][0]*tau_toty[bin_number]) - tau_fix[bin_number]) / (tau_fix[bin_number]*tau_fix_err_fac))**2
    
        else:
            p = 0
            for bin_number in num_bins:
                p -= 0.5 * ((params[model.config.par_slice(f"scale_nbsDstst_{bin_number}")][0] -constraints[f'scale_nbsDstst_{bin_number}'][0]) / (constraints[f'scale_nbsDstst_{bin_number}'][1])) ** 2
                p -= 0.5 * ((params[model.config.par_slice(f"scale_nbDstst_{bin_number}")][0] -constraints[f'scale_nbDstst_{bin_number}'][0]) / (constraints[f'scale_nbDstst_{bin_number}'][1])) ** 2
                p -= 0.5 * ((params[model.config.par_slice(f"scale_nDstDX_{bin_number}")][0] -constraints[f'scale_nDstDX_{bin_number}'][0]) / (constraints[f'scale_nDstDX_{bin_number}'][1])) ** 2
                mu_tau1 =  params[model.config.par_slice(f"scale_nDstTau_{bin_number}")][0] * full_templates[f"nDstTau_{bin_number}"]
                mu_tau2 =  params[model.config.par_slice(f"scale_nDstMu_{bin_number}")][0] * full_templates[f"nDstMu_{bin_number}"]
                p -= 0.5 * (((mu_tau1/mu_tau2) - constraints[f'MuTau_ratio_{bin_number}'][0]) / (constraints[f'MuTau_ratio_{bin_number}'][1]))**2

        return float(np.sum(logL) + p)
    init_pars = []
    for pp in model.config.par_order:
        if 'nbsDstst' in pp or 'nbDstst' in pp or 'nDstDX' in pp:
            init_pars.append(constraints[pp][0])
        else:
            init_pars.append(indep_res[pp][0])
    npars = len(init_pars)

    bounds = model.config.suggested_bounds() 
    init_pars = np.array(init_pars)


    def fcn(*pars):
        return -constrained_loglike(pars)
    name = model.config.par_order
 
    m = Minuit(fcn, *init_pars, name=name)
    for i in name:
        m.limits[i] = [0.0, upper_bounds[i]]
    m.errordef = 0.5 
    #m.tol = 0.001

    m.migrad()
    #m.hesse()
    param_dict = {}
    params = []
    params_ref = []
    for k in model.config.par_order:
        params_ref.append(m.values[k])
        param_dict[k] = [m.values[k], m.errors[k]]
        params.append(np.array([m.values[k], m.errors[k]]))
    params = np.array(params)
    dnll_0 = m.fval
    print(dnll_0)
    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.scipy_optimizer())


    limit_dict = {}


    limit_dict["bestfit_params"] = param_dict
    limit_dict["dnll_0"] = dnll_0
  

    limit_dict["fitted_yields"] = get_fitted_yields(wsp, params[:, 0])
    limit_dict["nominal_yields"] = {key : list(val.flatten()) for key, val in zip(wsp.model().config.samples, wsp.model().nominal_rates[0])}

    #limit_dict["min_nll"] = bestfit_params[w.model().config.poi_index]
    limit_dict["min_nll"] = params[:, 0]
    print('Param dict: ',param_dict)
    print('params ref: ', params_ref)


    return limit_dict
def constrain_ws_simultaneous(wsp, ratio, ratio_e, num_bins, toy, year, quad, linear):
    print('RATIO: ',ratio)
    print('RATIO e: ', ratio_e)
    model = wsp.model()
    data = wsp.data(model)
    wsmu_constraint = {}
    ddx_constraint = {}
    bs_toty = []
    b_toty = []
    sig_toty = []
    tau_toty = []
    yields_d = []
    yields_w = []
    yields_ddx = []
    
    for bin_number in range(num_bins):
        yw = 0
        yd = 0

        with open(fit_setup.template_paths['1D'] + f'{year}/' + f'data_WSMu_bin{bin_number}_raw.dat', "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                yw += float(v)


        with open(fit_setup.template_paths['1D'] + f'{year}/' + f'{fit_setup.data_path}{bin_number}_raw.dat', "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                v, e = line.split()
                yd +=float(v)
        wsmu_constraint[f'ratio_{bin_number}'] = yw/yd
        yields_w.append(yw)
        yields_d.append(yd)

    print(wsmu_constraint)
    print('len data ', len(data))
    def constrained_loglike(params):

        logL = model.logpdf(params,data)
        mu1 = 0
        mu2 = 0
        mu_tau1 = 0
        mu_tau2 = 0
        sum_bs = 0
        sum_tau = 0
        
        if toy:
            
            p = 0
            for bin_number in range(num_bins):
                mu1 += params[model.config.par_slice(f"scale_nWSMu_{bin_number}")][0] * yields_w[bin_number]
                mu2 += yields_w[bin_number]

            p -= 0.5 * (((mu1/mu2) - ratio) / ratio_e) ** 2

        else:


            p = 0
            for bin_number in range(num_bins):
                mu1 = params[model.config.par_slice(f"scale_nWSMu_{bin_number}")][0] * yields_w[bin_number]
                mu2 = yields_w[bin_number]

                p -= 0.5 * (((mu1/mu2) - ratio) / ratio_e) ** 2

        return float(np.sum(logL))
    init_pars = model.config.suggested_init()
    print('par_order: ', model.config.par_order)
    print('init pars: ',init_pars)
    npars = len(init_pars)

    bounds = model.config.suggested_bounds() 
    init_pars = np.array(init_pars)
    def fcn(*pars):
        return -constrained_loglike(pars)
    name = model.config.par_order
    m = Minuit(fcn, *init_pars, name=name)
    for i in name:
        m.limits[i] = [0.0, 100.0]
    m.errordef = 0.5 
    #m.tol = 0.001
    m.migrad()
    m.hesse()
    param_dict = {}
    params = []

    params_ref = []
    for k in model.config.par_order:
        params_ref.append(m.values[k])
        param_dict[k] = [m.values[k], m.errors[k]]
        params.append(np.array([m.values[k], m.errors[k]]))
    params = np.array(params)
    dnll_0 = m.fval
    print(dnll_0)
    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.scipy_optimizer())


    limit_dict = {}


    limit_dict["bestfit_params"] = param_dict
    limit_dict["dnll_0"] = dnll_0
  

    limit_dict["fitted_yields"] = get_fitted_yields(wsp, params[:, 0])
    limit_dict["nominal_yields"] = {key : list(val.flatten()) for key, val in zip(wsp.model().config.samples, wsp.model().nominal_rates[0])}
    limit_dict["min_nll"] = params[:, 0]
    print('Param dict: ',param_dict)
    print('params ref: ', params_ref)
    
    if quad or linear:
        if quad:
            smearlab = 'quad'
        else:
            smearlab = 'linear'
        def set_resolution(model, base_pars, value):
            pars_ = np.array(base_pars, copy=True)
            idx = model.config.par_order.index("Resolution")
            pars_[idx] = value
            return pars_
        fitted   = model.expected_actualdata(set_resolution(model, params_ref, limit_dict['bestfit_params']['Resolution'][0]))
        print(fitted)
        print(len(fitted))
        sample_contributions = {}

        for sample in model.config.samples:
            # print(sample)
            pars = params_ref.copy()
            for i, parname in enumerate(model.config.par_order):
                #print(parname)
                if parname not in ['scale_' + sample.replace('_w',''), 'Resolution']:
                        pars[i] = 0.0  # remove other samples

            pred = model.expected_actualdata(pars)
            sample_contributions[sample] = pred

        grouped = defaultdict(list)

        for name, arr in sample_contributions.items():
            idx = int(name.split("_")[-1])  # extract 0–4
            grouped[idx].append((name.replace(f'_w_{idx}',''), arr))
        n_per_group = 20

        for idx in sorted(grouped):
            yields = []
            errs = []
            bin_edges = fit_setup.mmiss2_bin_edges[f'1D_{fit_setup.n_w_bins}'][idx]
            bin_centers = (np.array(bin_edges[:-1]) + np.array(bin_edges[1:])) / 2
            bin_widths = np.diff(bin_edges)
            with open(fit_setup.template_paths['1D'] + f'{year}/' + f'{fit_setup.data_path}{idx}_raw.dat', "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    v, e = line.split()
                    yields.append(float(v))
                    errs.append(float(e))
            data= np.array(yields)
            data_err = np.array(errs)
            start = idx * n_per_group
            end = (idx + 1) * n_per_group

            bins = np.arange(n_per_group)
            bottom = np.zeros(n_per_group)
            #full_err = np.zeros(n_per_group)
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True ,gridspec_kw={'height_ratios': [3, 1]})

            for name, arr in reversed(grouped[idx]):
                vals = arr[start:end]
                ax1.bar(bin_centers, vals, width=bin_widths,bottom=bottom, color = fit_setup.component_colors[name],label=fit_setup.component_labels[name]+F': {np.sum(vals):.0f}')
                bottom += vals
            ax1.errorbar(bin_centers, data, yerr = data_err, color = 'black', label = 'data' + f': {np.sum(data):.0f}', fmt = 'o',  ms = 2, capsize = 1)
            ax1.set_ylim(0,1.1*bottom.max())
            ax1.set_xlabel("Bin")
            ax1.set_ylabel("Events")
            ax1.legend(reverse = True)
            
            pulls = (bottom- data)/np.sqrt( bottom + data_err**2)
            #pulls = (bot- data)/np.sqrt(data_err**2)

            ax2.axhline(0, color="black", linestyle="--", linewidth=0.8)
            ax2.errorbar(bin_centers, pulls, yerr=np.ones_like(pulls), fmt="o", color="black")
            ax2.set_ylim(-5,5)
            ax2.set_ylabel("Pulls")
            plt.savefig(f'check_{smearlab}_bin{idx}_{year}.pdf')
            plt.clf()

    return limit_dict

def get_fit_params_errors(w, bestfit_params=None):
    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.minuit_optimizer(tolerance=1e-3))
    params = pyhf.infer.mle.fit(
        w.data(w.model()),
        w.model(),
        return_uncertainties=True,
        return_correlations=True,
        init_pars=bestfit_params,
    )
    pyhf.set_backend(pyhf.tensorlib, pyhf.optimize.scipy_optimizer())
    return params



def get_fitted_yields(w, params):
    #taken from _MainModel.expected_data at
    #https://github.com/scikit-hep/pyhf/blob/main/src/pyhf/pdf.py
    model = w.model()
    deltas, factors = model._modifications(params)
    allsum = np.concatenate(deltas + [model.nominal_rates])
    nom_plus_delta = np.sum(allsum, axis=0)
    nom_plus_delta = np.reshape(nom_plus_delta, (1,) + np.shape(nom_plus_delta))
    allfac = np.concatenate(factors + [nom_plus_delta])
    newbysample = np.product(allfac, axis=0)
    newresults = np.sum(newbysample, axis=0)
    return {key : list(val.flatten()) for key, val in zip(model.config.samples, newbysample)}



def plot_fitted_dist(bin_idx, year, fit_params, extra_label, dir,  dim = '1D', fix_dsttau = False, samples = ["nbDstst", "nbsDstst", "nDstst_higher", "nDstTau", "nDstDX", 'nWSMu', 'nWSPi', "nDstMu"],jsons_dir = fit_setup.template_path, toy = False):
    label = f"{year}_w_{bin_idx}{extra_label}"
    out_lab = label
    leg_lab = fit_setup.component_labels
    yields = []
    errs= []
    bin_edges = fit_setup.mmiss2_bin_edges[f'{dim}_{fit_setup.n_w_bins}'][bin_idx]

    bin_centers = (np.array(bin_edges[:-1]) + np.array(bin_edges[1:])) / 2
    bin_widths = np.diff(bin_edges)
    with open(fit_setup.template_paths[dim] + f'{year}/' + f'{fit_setup.data_path}{bin_idx}_raw.dat', "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            v, e = line.split()
            yields.append(float(v))
            errs.append(float(e))
    data= np.array(yields)
    data_err = np.array(errs)
    scale_data = np.sum(data)
    if fix_dsttau:
        out_lab += '_DstTau_fixed'
    if toy:
        templates = json.load(open(jsons_dir / f"toys/histos_for_fit_{out_lab}.json"))
        numbers = ''
        for char in reversed(out_lab):
            if char.isdigit():
                numbers = char + numbers
            else:
                break
        fracs = json.load(open(jsons_dir / f"toys/fractions_{numbers}.json"))
    else:
        templates = {}
        for k,v in fit_setup.component_templates.items():
            yields = []
            errs = []
            with open(fit_setup.template_paths[dim] + f'{year}/' + f'{v}{bin_idx}_raw.dat', "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    v, e = line.split()
                    yields.append(float(v))
                    errs.append(float(e))
                templates[k] =  np.array(yields)

    errs = {}
    colors = fit_setup.component_colors

    bot = np.zeros(len(bin_centers))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True ,gridspec_kw={'height_ratios': [3, 1]})
    for s in samples:

        if s == 'nDstTau' and fix_dsttau:
            fit_param = 1
            errs['nDstTau'] = np.array(templates['nDstTau_stat'])
        else:
            fit_param = fit_params[f'scale_{s}_{bin_idx}'][0]
            err = fit_params[f'scale_{s}_{bin_idx}'][1]
            errs[s] = fit_params[f'scale_{s}_{bin_idx}'][1]*np.array(templates[s])

        if toy:
            ax1.bar(bin_centers, fit_param*np.array(templates[s]), width = bin_widths, bottom=bot, color = colors[s], label=f'{leg_lab[s]}')#: {'%.3f'%(fit_param)}$\pm${'%.3f'%(err)} ({fracs[s]})')
        else:
            ax1.bar(bin_centers, fit_param*np.array(templates[s]), width = bin_widths, bottom=bot, color = colors[s], label=f'{leg_lab[s]}' + f': {(fit_param*np.sum(templates[s])):.0f}')#: {'%.0f'%(fit_param*np.sum(templates[s]))}')
        bot += fit_param*np.array(templates[s])
    #ax1.set_yscale('log')
    ax1.errorbar(bin_centers, data,xerr = 0.5*bin_widths, yerr= data_err, fmt = 'o', ms = 2, capsize = 1, color = 'black', label = f'data' + f': {(np.sum(data)):.0f}')#: {'%.0f'%(np.sum(templates['data']))}')
    ax1.legend(reverse = 'True', fontsize = 'small', frameon = False, loc = 'upper left')
    full_err = np.zeros(len(bin_centers))
    for s in samples:
        full_err += errs[s]**2

    pulls = (bot- data)/np.sqrt(data_err**2)

    ax2.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax2.errorbar(bin_centers, pulls, yerr=np.ones_like(pulls), fmt="o", color="black")
    ax2.set_ylim(-5,5)
    ax2.set_ylabel("Pulls")
    if dim == '5D':
        ax1.set_title(fit_setup.bin_labels_5D[bin_idx])
    else:
        ax1.set_title(f'BIN {bin_idx}')
    plt.savefig(dir / f'fitted_distributions_{out_lab}_with_WSMu_witherr.pdf')



def plot_fitted_dist_raw(bin_idx, year, fit_params, extra_label, dir,smearing,  dim = '1D', fix_dsttau = False, samples = ["nbDstst", "nbsDstst", "nDstst_higher", "nDstTau", "nDstDX", 'nWSMu', 'nWSPi', "nDstMu"],jsons_dir = fit_setup.template_path, toy = False, smear = False):
    label = f"{year}_w_{bin_idx}{extra_label}"
    out_lab = label
    leg_lab = fit_setup.component_labels
    yields = []
    errs= []
    bin_edges = fit_setup.mmiss2_bin_edges[f'{dim}_{fit_setup.n_w_bins}'][bin_idx]

    bin_centers = (np.array(bin_edges[:-1]) + np.array(bin_edges[1:])) / 2
    bin_widths = np.diff(bin_edges)
    with open(fit_setup.template_paths[dim] + f'{year}/' + f'{fit_setup.data_path}{bin_idx}_raw.dat', "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            v, e = line.split()
            yields.append(float(v))
            errs.append(float(e))
    data= np.array(yields)
    data_err = np.array(errs)
    scale_data = np.sum(data)
    if fix_dsttau:
        out_lab += '_DstTau_fixed'
    else:
        templates = {}
        temp_lab = ''
        for k,v in fit_setup.component_templates.items():
            if not smear and k != 'nWSMu':
                temp_lab = '_final'
            else:
                temp_lab = ''
            yields = []
            errs = []
            with open(fit_setup.template_paths[dim] + f'{year}/' + f'{v}{bin_idx}{temp_lab}_raw.dat', "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    v, e = line.split()
                    yields.append(float(v))
                    errs.append(float(e))
                templates[k] = np.array(yields)#  * np.diff(bin_edges)

    errs = {}
    colors = fit_setup.component_colors

    bot = np.zeros(len(bin_centers))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True ,gridspec_kw={'height_ratios': [3, 1]})
    for s in samples:

        if s == 'nDstTau' and fix_dsttau:
            fit_param = 1
            errs['nDstTau'] = np.array(templates['nDstTau_stat'])
        else:
            fit_param = fit_params['bestfit_params'][f'scale_{s}_{bin_idx}'][0]
            err = fit_params['bestfit_params'][f'scale_{s}_{bin_idx}'][1]*np.array(templates[s])
            errs[s] = fit_params['bestfit_params'][f'scale_{s}_{bin_idx}'][1]*np.array(templates[s])

        if toy:
            ax1.bar(bin_centers, fit_param*np.array(templates[s]), width = bin_widths, bottom=bot, color = colors[s], label=f'{leg_lab[s]}')#: {'%.3f'%(fit_param)}$\pm${'%.3f'%(err)} ({fracs[s]})')
        else:
            ax1.bar(bin_centers, fit_param*np.array(templates[s]), width = bin_widths, bottom=bot, color = colors[s], label=f'{leg_lab[s]}' + f': {(fit_param*np.sum(templates[s])):.0f}')#: {'%.0f'%(fit_param*np.sum(templates[s]))}')
        bot += fit_param*np.array(templates[s])
    #ax1.set_yscale('log')
    ax1.errorbar(bin_centers, data,xerr = 0.5*bin_widths, yerr= data_err, fmt = 'o', ms = 2, capsize = 1, color = 'black', label = f'data' + f': {(np.sum(data)):.0f}')#: {'%.0f'%(np.sum(templates['data']))}')
    ax1.legend(reverse = 'True', fontsize = 'small', frameon = False, loc = 'upper left')
    full_err = np.zeros(len(bin_centers))
    for s in samples:
        full_err += errs[s]**2
    for i, err in zip(range(len(bin_centers)), full_err):
        x = bin_edges[i]
        height = bot[i]
        ax1.add_patch(plt.Rectangle((x, height-np.sqrt(err)/2), bin_widths[i], np.sqrt(err), fill=True, facecolor = 'darkgrey', alpha = 0.4,
                                edgecolor='black', hatch='////', lw = 0))
    pulls = (bot- data)/np.sqrt( full_err + data_err**2)


    ax2.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax2.errorbar(bin_centers, pulls, yerr=np.ones_like(pulls), fmt="o", color="black")
    ax2.set_ylim(-5,5)
    ax2.set_ylabel("Pulls")
    if dim == '5D':
        ax1.set_title(fit_setup.bin_labels_5D[bin_idx])
    else:
        ax1.set_title(f'BIN {bin_idx}')
    plt.savefig(dir / f'fitted_distributions_{out_lab}_with_WSMu_witherr_raw.pdf')




def plot_fitted_dist_raw_control(bin_idx, year, fit_params, extra_label, dir,smearing,  dim = '1D', fix_dsttau = False, samples = ["nbDstst", "nbsDstst", "nDstst_higher", "nDstTau", "nDstDX", 'nWSMu', 'nWSPi', "nDstMu"],jsons_dir = fit_setup.template_path, toy = False, smear = False):
    label = f"{year}_w_{bin_idx}{extra_label}"
    out_lab = label
    leg_lab = fit_setup.component_labels
    yields = []
    errs= []
    bin_edges = fit_setup.mmiss2_bin_edges[f'{dim}_{fit_setup.n_w_bins}'][bin_idx]

    bin_centers = (np.array(bin_edges[:-1]) + np.array(bin_edges[1:])) / 2
    bin_widths = np.diff(bin_edges)
    with open(fit_setup.template_paths[dim] + f'{year}/' + f'{fit_setup.data_path}{bin_idx}_raw.dat', "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            v, e = line.split()
            yields.append(float(v))
            errs.append(float(e))
    data= np.array(yields)
    data_err = np.array(errs)
    scale_data = np.sum(data)
    if fix_dsttau:
        out_lab += '_DstTau_fixed'
    else:
        templates = {}
        temp_lab = ''
        for k,v in fit_setup.component_templates.items():
            if not smear and k != 'nWSMu':
                temp_lab = ''
            else:
                temp_lab = ''
            yields = []
            errs = []
            with open(fit_setup.template_paths[dim] + f'{year}/' + f'{v}{bin_idx}{temp_lab}_raw.dat', "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    v, e = line.split()
                    yields.append(float(v))
                    errs.append(float(e))
                templates[k] = np.array(yields)#  * np.diff(bin_edges)
    errs = {}
    colors = fit_setup.component_colors
    templates['dstst'] = templates['nbsDstst'] + templates['nbDstst'] + templates['nDstst_higher']
    bot = np.zeros(len(bin_centers))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True ,gridspec_kw={'height_ratios': [3, 1]})
    for s in samples:

        if s == 'nDstTau' and fix_dsttau:
            fit_param = 1
            errs['nDstTau'] = np.array(templates['nDstTau_stat'])
        else:
            fit_param = fit_params[f'scale_{s}_{bin_idx}'][0]
            err = fit_params[f'scale_{s}_{bin_idx}'][1]*np.array(templates[s])
            errs[s] = fit_params[f'scale_{s}_{bin_idx}'][1]*np.array(templates[s])

        if toy:
            ax1.bar(bin_centers, fit_param*np.array(templates[s]), width = bin_widths, bottom=bot, color = colors[s], label=f'{leg_lab[s]}')#: {'%.3f'%(fit_param)}$\pm${'%.3f'%(err)} ({fracs[s]})')
        else:
            if 'Dstst' in s:
                if s == 'nDstst_higher':
                    dststy = fit_param*np.sum(templates['dstst'])
                    ax1.bar(bin_centers, fit_param*np.array(templates[s]), width = bin_widths, bottom=bot, color = 'violet', label='$B \\to D^{**}\\mu\\nu$' + f': {(dststy):.0f}')
                else:
                    ax1.bar(bin_centers, fit_param*np.array(templates[s]), width = bin_widths, bottom=bot, color = 'violet')
            elif s == 'nDstTau':
                ax1.bar(bin_centers, fit_param*np.array(templates[s]), width = bin_widths, bottom=bot, color = colors[s])#: {'%.0f'%(fit_param*np.sum(templates[s]))}')
            else:   
                ax1.bar(bin_centers, fit_param*np.array(templates[s]), width = bin_widths, bottom=bot, color = colors[s], label=f'{leg_lab[s]}' + f': {(fit_param*np.sum(templates[s])):.0f}')#: {'%.0f'%(fit_param*np.sum(templates[s]))}') 
        bot += fit_param*np.array(templates[s])
    #ax1.set_yscale('log')
    ax1.errorbar(bin_centers, data,xerr = 0.5*bin_widths, yerr= data_err, fmt = 'o', ms = 2, capsize = 1, color = 'black', label = f'data' + f': {(np.sum(data)):.0f}')#: {'%.0f'%(np.sum(templates['data']))}')
    ax1.legend(reverse = 'True', fontsize = 'small', frameon = False, loc = 'upper left')
    full_err = np.zeros(len(bin_centers))
    for s in samples:
        full_err += errs[s]**2
    for i, err in zip(range(len(bin_centers)), full_err):
        x = bin_edges[i]
        height = bot[i]
        ax1.add_patch(plt.Rectangle((x, height-np.sqrt(err)/2), bin_widths[i], np.sqrt(err), fill=True, facecolor = 'darkgrey', alpha = 0.4,
                                edgecolor='black', hatch='////', lw = 0))
    pulls = (bot- data)/np.sqrt( full_err + data_err**2)
    #pulls = (bot- data)/np.sqrt(data_err**2)

    ax2.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax2.errorbar(bin_centers, pulls, yerr=np.ones_like(pulls), fmt="o", color="black")
    ax2.set_ylim(-5,5)
    ax2.set_ylabel("Pulls")

    if dim == '5D':
        ax1.set_title(fit_setup.bin_labels_5D[bin_idx])
    else:
        ax1.set_title(f'BIN {bin_idx}')
    plt.savefig(dir / f'fitted_distributions_{out_lab}_with_WSMu_witherr_raw.pdf')

