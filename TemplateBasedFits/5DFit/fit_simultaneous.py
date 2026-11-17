#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to make the pyhf json to run fits
"""

import numpy as np
import pyhf
#import copy
import json
from pathlib import Path
import matplotlib.pyplot as plt
import fit_meth
import fit_setup

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-y",
        "--year",
        help="""Which year dataset to use""",
        choices=['2016', '2017', '2018', 'Run2'],
        type=str,
        default='2016',
    )
    parser.add_argument(
        "-d",
        "--dim",
        help="""1D or 5D fit""",
        choices=['1D', '5D'],
        type=str,
        default='1D',
    )
    parser.add_argument(
        "-i",
        "--include",
        help="""Which samples to use in toys""",
        action='append',
        default= [],
    )

    parser.add_argument(
        "--eff_scale", help="scale the signal efficiency", type=float, default=1
    )
    parser.add_argument(
        "--mu_scale",
        help="scale the signal in the expected data",
        type=float,
        default=1,
    )
    parser.add_argument(
        "--fix_tau", help="fix tau contribution", action="store_true"
    )
    parser.add_argument(
        "--no_syst",
        help="Do not consider systematic uncertainties",
        action="store_false",
        default = True
    )
    parser.add_argument(
        "--max",
        help="upper limit of the mu range scanned (default=10)",
    )
    parser.add_argument(
        "--min",
        help="upper limit of the mu range scanned (default=0)",
    )
    parser.add_argument(
        "-v",
        "--values",
        help="Give json file with the pyhf workspace to compute the limit",
    )
    parser.add_argument(
        "-n", help="Points to compute the profile likelihood", type=int, default=20
    )
    parser.add_argument(
        "--toy", help="use Poisson-fluctuated data",  action="store_true"
    )
    parser.add_argument("--no_weights", help="Use no weights", action="store_false", default = True)
    parser.add_argument("--constrain_feeddown", help="Use no weights", action="store_true")
    parser.add_argument("--smear", help="whether to float smearing", action="store_true")

    args = parser.parse_args()
    dim = args.dim
    quad = False
    linear = False
    constrain_feeddown = args.constrain_feeddown
    year = args.year
    extra_label = dim
    smear = args.smear
    ii= args.include
    if ii == []:
        include = ["nbDstst", "nbsDstst", "nDstst_higher", "nDstTau", "nDstDX", 'nWSMu']
    else:
        include = ii
    if not args.no_weights:
        extra_label += f"_Weights"
    toy = args.toy
    fix_tau = args.fix_tau
  
    bin_idx = 0
    label = f"{year}_{bin_idx}{extra_label}"
    cl = 0.9
    out_label = label

    if args.mu_scale != 1:
        out_label += f"_mu{args.mu_scale}"
        
    this_dir = Path(__file__).parent.resolve()
    out_dir = this_dir / f"fit_results_{dim}"
    if bool(toy):
        toy_str = ''
        for k,v in fit_setup.fit_fixed_fracs.items():
            aa = k.replace('scale_','')
            toy_str += aa + '_' + str(v) + '_'
            toy_str = toy_str[:-1]
        figdir = out_dir / f"results_toys_{toy_str}"
    else:
        figdir = out_dir / f"fit_plots_{dim}"
    specdir = out_dir / f"spec"
    figdir.mkdir(parents=True, exist_ok=True)
    specdir.mkdir(parents=True, exist_ok=True)
    if dim == '5D':
        nbins = fit_setup.nbins_4D['w_Pred'] * fit_setup.nbins_4D['costhetal_reco'] * fit_setup.nbins_4D['costhetad_reco'] * fit_setup.nbins_4D['chi_reco']
        
    else:
        nbins = fit_setup.n_w_bins
    print(label)
    print(include)
    
    if args.values:
        spec = json.load(open(args.values))

    else: 
        spec = {
            "channels": [],
            "observations": [],
            "measurements": [],
            "version": "1.0.0",
        }

        #spec_nosysts = copy.deepcopy(spec)
        poi_name = 'mu_sig'
        fit_meth.add_channel_q2_bins_simultaneous(
            spec,
            year,
            wbin = bin_idx,
            dim = dim,
            extra_label=extra_label,
            #lumi=args.lumi,
            fix_dsttau= fix_tau,
            mu_scale=args.mu_scale,
            eff_scale=args.eff_scale,
            syst=not args.no_syst,
            toy = bool(toy),
            incl = include,
            smear = smear,
            quad = quad,
            linear = linear
 
        )

        fit_meth.add_measurement_q2_bins_simultaneous(spec, nbins = nbins,max_range=args.max, incl = include)
    print(spec)
    w = pyhf.workspace.Workspace(spec)
    if args.min or args.max:
        mu_tests = np.linspace(
            float(args.min) if args.min else 0,
            float(args.max) if args.max else 10,
            args.n,
        )
    else:
        if constrain_feeddown:
            sig_mu = 1
            for k in include: 
                sig_mu -= fit_setup.fit_fixed_fracs[k]
            ratio_tau = fit_setup.fit_fixed_fracs['nDstTau']/sig_mu
            print('nDstMu', sig_mu)
            limit_dict = fit_meth.constrain_bs_simultaneous(w, nbins, toy, year, smear = smear)            
            fit_params = limit_dict["bestfit_params"]

        else:

            limit_dict = fit_meth.get_fitted_params(w)
            fit_params = limit_dict["bestfit_params"]

    if not toy:
        for bin_idx in range(nbins):
            fit_meth.plot_fitted_dist_raw(bin_idx, year, fit_params, extra_label, figdir, smearing= smear, dim = dim,toy = bool(toy), samples = include + ['nDstMu'], smear = smear)
        f = open(figdir / f'fit_result_{extra_label}_{year}.txt', 'w')
        samples = ['nDstMu'] + include
        for i in fit_params.keys():
            fit_param = fit_params[i][0]
            err = fit_params[i][1]
            a = str(i) + " "
            a += str(fit_param)
            a += " "
            a += str(err)
            f.write(a+"\n")
        f.close()
    else:

        random_seed = np.random.RandomState()
        toy_id = np.random.randint(0,1e8)
        f = open(figdir / f'toy_result_{toy_id}.txt', 'w')
        #samples = ["nbDstst", "nbsDstst", "nDstst_higher", "nDstTau", "nDstDX", 'nWSMu', 'nWSPi', "nDstMu"]
        samples = ['nDstMu'] + include
        for i in fit_params.keys():
            fit_param = fit_params[i][0]
            err = fit_params[i][1]
            a = str(i) + " "
            a += str(fit_param)
            a += " "
            a += str(err)
            f.write(a+"\n")
        f.close()

        
        f = open (figdir / f'toy_list_bin{bin_idx}.txt','a')
        f.write(str(toy_id) + "\n")
        f.close()
