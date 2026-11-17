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
        "-b",
        "--binIdx",
        help="""Which bin should the m_corr fit be performed in""",
        choices=list(range(225)) + [9999],
        type=int,
        default=1,
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
        choices=['1D', '5D', 'control'],
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
        "--fix_tau", help="fix the tau contirbution", action="store_true"
    )
    parser.add_argument(
        "--no_syst",
        action="store_false",
        default = True
    )

    parser.add_argument("--cls", help="Compute also CLs", action="store_true")
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
    smear = args.smear
    constrain_feeddown = args.constrain_feeddown
    year = args.year
    extra_label = dim
    ii= args.include
    if ii == []:
        include = ["nbDstst", "nbsDstst", "nDstst_higher", "nDstTau", "nDstDX", 'nWSMu']#, 'nWSPi']
    else:
        include = ii
    if not args.no_weights:
        extra_label += "_Weights"
    toy = args.toy
    fix_tau = args.fix_tau
  
    bin_idx = args.binIdx
    label = f"{year}_{bin_idx}{extra_label}"
    cl = 0.9
    out_label = label

    if args.mu_scale != 1:
        out_label += f"_mu{args.mu_scale}"
        
    this_dir = Path(__file__).parent.resolve()
    out_dir = this_dir / f"fit_results_{dim}"
    if bool(toy):
        figdir = out_dir / f"results_toys_5D"
    else:
        figdir = out_dir / f"fit_plots_{dim}"
    specdir = out_dir / f"spec"
    figdir.mkdir(parents=True, exist_ok=True)
    specdir.mkdir(parents=True, exist_ok=True)
   
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
        fit_meth.add_channel_q2_bins(
            spec,
            bin_idx,
            year,
            dim = dim,
            extra_label=extra_label,
            #lumi=args.lumi,
            fix_dsttau= fix_tau,
            mu_scale=args.mu_scale,
            eff_scale=args.eff_scale,
            syst=not args.no_syst,
            toy = bool(toy),
            incl = include,
            smear = smear
 
        )
            
        fit_meth.add_measurement_q2_bins(spec, bin_idx, fix_dsttau= fix_tau, dim = dim,max_range=args.max, incl = include)
    w = pyhf.workspace.Workspace(spec)

    if args.min or args.max:
        mu_tests = np.linspace(
            float(args.min) if args.min else 0,
            float(args.max) if args.max else 10,
            args.n,
        )
    
    else:
        if constrain_feeddown:
            eff_bsdstst = 1
            eff_bdstst = 1
            wbin = bin_idx // 45
            mu_tests = fit_meth.constrain_5d_indep(w, bin_idx = bin_idx, wbin = wbin, toy = toy, year = year)

        else:
            mu_tests = fit_meth.auto_mu_tests(w, bin_idx, args.n)
            limit_dict = fit_meth.get_limit_dict(
        w, mu_tests=mu_tests, cl=cl, n_points=args.n, no_systs=args.s
    )
            fit_params = {
        name: val
        for name, val in zip(w.model().config.par_names, limit_dict["bestfit_params"])
    }
            del fit_params[w.model().config.par_names[w.model().config.poi_index]]
            fit_params_errors = {
                name: val[1]
                for name, val in zip(
                    w.model().config.par_names, limit_dict["bestfit_params_errors"]
                )
            }

        #print(mu_tests)

    if not toy:

        fit_meth.plot_fitted_dist_raw(bin_idx, year, mu_tests, extra_label, figdir, dim = dim,toy = bool(toy), samples = include + ['nDstMu'], smearing = False)
        f = open(figdir / f'fit_result_{extra_label}_{bin_idx}_{year}.txt', 'w')
        samples = ['nDstMu'] + include
        for i in mu_tests['bestfit_params'].keys():
            fit_param = mu_tests['bestfit_params'][i][0]
            err = mu_tests['bestfit_params'][i][1]
            a = str(i) + " "
            a += str(fit_param)
            a += " "
            a += str(err)
            f.write(a+"\n")
        f.close()
    else:
       
        f = open (figdir / f'toy_list_bin{bin_idx}.txt','a')
        f.write(str(mu_tests[0]) + " " + str(mu_tests[1]) + " " + str(mu_tests[2])   + "\n")
        #f.write(str(toy_id) + "\n")
        f.close()
        np.savez(figdir / f"cov_{bin_idx}.npz",
         covariance=mu_tests[3],
         parameters=mu_tests[4])
