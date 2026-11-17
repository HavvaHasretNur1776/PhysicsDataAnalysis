#!/usr/bin/env python3
# -*- coding: utf-8 -


import json
import os
import copy
import sys
import numpy as np
import ROOT as r
import fit_setup
import pandas as pd
import matplotlib.pyplot as plt





n_w_bins = 5
year = str(sys.argv[1])
dim = str(sys.argv[2])
control = str(sys.argv[3]) == 'control'
mmiss = '_final'
if control:
    edges = 'control'
else:
    edges = dim + f'_{n_w_bins}'
template_path = fit_setup.template_paths[dim]
postfit_bins = 20
postfit_edges = {
    #'mu_PE' : np.linspace(0,100000, postfit_bins),
    #'m_corr_reco' : np.linspace(2000,9000, postfit_bins),
    # 'b0_P_Pred' : np.linspace(0, 400000, postfit_bins),
    #'q2GeV_reco_regreg' : np.linspace(0, 11, 5),
}
json_dir = fit_setup.jsons_dir
bin_vars= fit_setup.angular_vars
bin_file = fit_setup.bin_file
nbins = fit_setup.nbins_4D
ranges = fit_setup.ranges_4D
bins_w = fit_setup.w_bins_4D
if dim == '1D':
    bins = fit_setup.mmiss2_bin_edges[edges]
elif dim == '5D':
    bins = json.load(open(json_dir + bin_file))

v0 = bin_vars[0]
v1 = bin_vars[1]
v2 = bin_vars[2]
v3 = bin_vars[3]

WSMu_collapsed_vars = ['ctd_reco_regreg', 'chi_reco_regreg']

if dim == '5D':
    bin_labels = {}
    binno = 0
    for bv0 in range(len(bins[v0])-1):
        for bv1 in range(len(bins[v1])-1):
            for bv2 in range(len(bins[v2])-1):
                for bv3 in range(len(bins[v3])-1):
                    bin_labels[binno] = f'{bv0}_{bv1}'
                    binno+=1
else:
    binno = len(fit_setup.w_edges[n_w_bins]) -1

tuples = {
  'MC' : {
    'nDstMu': 'b0_dstmunu',
    'nDstTau': 'b0_dsttaunu',
    'nDstDX' : {
        "b+_dstdx": 'b+_dstdx',
        "b0_dstdx": 'b0_dstdx',
    },
    'nbDstst' : {
        "b0_dststmu": 'b0_dststmu',
        "b0_dststmu_415": 'b0_dststmu',
        "b0_dststmu_10413": 'b0_dststmu',
        "b0_dststmu_20413": 'b0_dststmu',

        "b+_dststmu": 'b+_dststmu',
        "b+_dststmu_425": 'b+_dststmu',
        "b+_dststmu_10423": 'b+_dststmu',
        "b+_dststmu_20423": 'b+_dststmu',
    },
    'nbsDstst' : {
        "bs_dststmu_10433": "bs_dststmu",
        "bs_dststmu_435": "bs_dststmu",

    },
    'nDstst_higher' : {
        "b+_dststmu_heavy": "b+_dststmu_heavy",
        "b0_dststmu_heavy": "b0_dststmu_heavy",

    },
  },
'data' : {
    'nWSMu': 'data_WSMu',
    'nWSPi': 'data_WSPi',# 'WSMuWSPi_BestVariables_angles.root'],
    'data' :  "data_RS"}, }





template_label = {
'nbsDstst': 'bs_Dstst',
'nDstTau': 'b0_DstTau',
'nDstDX': 'DstDx',
'nbDstst': 'b_Dstst',
'nDstst_higher': 'Dstst_higher',
'nDstMu': 'signal',
'nWSMu': 'data_WSMu',
'nWSPi': 'data_WSPi',
'data' :  "data",
 }

mc_weight_base = [
    "PIDCalibEff",
    "w_trk_mu",
    "w_trk_k",
    "w_trk_pi",
    "w_trk_spi",
    "w_l0_tosORtis_emu",
    "w_spdcorr",
    "w_bkin_tot",
]

weights_hammer = {
    'nDstMu': "SchemeBGL",
    'nDstTau': "SchemeBGL",
 
    'nbDstst' : {
        "b0_dststmu_415": "SchemeBLR",
        "b0_dststmu_10413": "SchemeBLR",
        "b0_dststmu_20413": "SchemeBLR",
        
        "b+_dststmu_425": "SchemeBLR",
        "b+_dststmu_10423": "SchemeBLR",
        "b+_dststmu_20423": "SchemeBLR",
    },
    'nbsDstst' : {
        "bs_dststmu_10433": "SchemeBLR",
        "bs_dststmu_435": "SchemeBLR",

    },
}

tm_strings = {
    'nDstMu': 'abs(b0_TRUEID) == 511 and abs(mu_MC_MOTHER_ID) == 511 and abs(dst_TRUEID) == 413',
    'nDstTau': 'abs(b0_TRUEID) == 511 and abs(mu_MC_MOTHER_ID) == 15 and abs(dst_TRUEID) == 413',
    'nDstDX' : {
        "b+_dstdx": 'abs(b0_TRUEID) == 521 and abs(mu_MC_MOTHER_ID) == 431 and abs(dst_TRUEID) == 413',
        "b0_dstdx": 'abs(b0_TRUEID) == 511 and abs(mu_MC_MOTHER_ID) == 431 and abs(dst_TRUEID) == 413',
    },
    'nbDstst' : {
        "b0_dststmu_415": 'abs(b0_TRUEID) == 511 and abs(mu_MC_MOTHER_ID) == 511 and abs(dst_TRUEID) == 413 and abs(b0_TrueHadron_D0_ID) == 415',
        "b0_dststmu_10413": 'abs(b0_TRUEID) == 511 and abs(mu_MC_MOTHER_ID) == 511 and abs(dst_TRUEID) == 413 and abs(b0_TrueHadron_D0_ID) == 10413',
        "b0_dststmu_20413": 'abs(b0_TRUEID) == 511 and abs(mu_MC_MOTHER_ID) == 511 and abs(dst_TRUEID) == 413 and abs(b0_TrueHadron_D0_ID) == 20413',
        
        "b+_dststmu_425": 'abs(b0_TRUEID) == 521 and abs(mu_MC_MOTHER_ID) == 521 and abs(dst_TRUEID) == 413 and abs(b0_TrueHadron_D0_ID) == 425',
        "b+_dststmu_10423": 'abs(b0_TRUEID) == 521 and abs(mu_MC_MOTHER_ID) == 521 and abs(dst_TRUEID) == 413 and abs(b0_TrueHadron_D0_ID) == 10423',
        "b+_dststmu_20423": 'abs(b0_TRUEID) == 521 and abs(mu_MC_MOTHER_ID) == 521 and abs(dst_TRUEID) == 413 and abs(b0_TrueHadron_D0_ID) == 20423',
    },
    'nbsDstst' : {
        "bs_dststmu_10433": 'abs(b0_TRUEID) == 531 and abs(mu_MC_MOTHER_ID) == 531 and abs(dst_TRUEID) == 413 and abs(b0_TrueHadron_D0_ID) == 10433',
        "bs_dststmu_435": 'abs(b0_TRUEID) == 531 and abs(mu_MC_MOTHER_ID) == 531 and abs(dst_TRUEID) == 413 and abs(b0_TrueHadron_D0_ID) == 435',

    },
    'nDstst_higher' : {
        "b+_dststmu_heavy": 'abs(b0_TRUEID) == 521 and abs(mu_MC_MOTHER_ID) == 521 and abs(dst_TRUEID) == 413 and abs(b0_TrueHadron_D0_ID) == 425',
        "b0_dststmu_heavy": 'abs(b0_TRUEID) == 511 and abs(mu_MC_MOTHER_ID) == 511 and abs(dst_TRUEID) == 413 and abs(b0_TrueHadron_D0_ID) == 415',
    },
}

rel_BFs = {
 
    'nDstDX' : {
        "b+_dstdx": [0.5, 0],
        "b0_dstdx": [0.5, 0],
    },
    'nbDstst' : {
        "b0_dststmu_415": [0.5 * 6.6 * 1e-4, 0.5 * 1.1 * 1e-4],#0.5 factors for D**- -> D*- isospin symmetry, keep original measurement errors
        "b0_dststmu_10413": [0.5 * 2.85 * 1e-3, 0.5 * 0.25*1e-3],
        "b0_dststmu_20413": [0.5 * 2.5 *1e-3, 0.5 * 0.5 * 1e-3],
        
        "b+_dststmu_425": [1.06*1e-3, 0.28*1e-3],
        "b+_dststmu_10423": [2.84 * 1e-3, 0.17 * 1e-3],
        "b+_dststmu_20423": [1.7 * 1e-3, 0.6 * 1e-3],
    },
    'nbsDstst' : {
        "bs_dststmu_10433": [2 * 2.7 * 1e-3, 2 * 0.7 * 1e-3],
        "bs_dststmu_435": [4.3 * 1e-5, 2.24 * 1e-5],

    },
    'nDstst_higher' : {
        "b+_dststmu_heavy": [0.5 * 6.6 * 1e-4, 0.5 * 1.1 * 1e-4],
        "b0_dststmu_heavy": [1.06 * 1e-3, 0.18 * 1e-3],
    },
}




component_labels = {

    "n415": r"$B^0 \to D_2^*(2460)^+ \mu^- \nu_\mu$",
    "n10423": r"$B^+ \to D_1(2420)^0 \mu^- \nu_\mu$",
    "n20413": r"$B^0 \to D_1(H)^+ \mu^- \nu_\mu$",
    "n435": r"$B_s \to D_{s2}^*(2573)^+ \mu^- \nu_\mu$",
    "n425": r"$B^+ \to D_2^*(2460)^0 \mu^- \nu_\mu$",
    "nBdDststHigher": r"$B^0 \to D^{**}(Higher) \mu^- \nu_\mu$",
    "nBuDststHigher": r"$B^+ \to D^{**}(Higher) \mu^- \nu_\mu$",
    "n10433": r"$B_s \to D_{s1}(2536)^+ \mu^- \nu_\mu$",
    "n10413": r"$B^0 \to D_1(2420)^+ \mu^- \nu_\mu$",
    "n20423": r"$B^+ \to D_1(2430)^0 \mu^- \nu_\mu$",
}




dec_fracs = {
    'nDstDX' : {
        "b+_dstdx": 1,
        "b0_dstdx": 1,
    },
    'nbDstst' : {
        "b0_dststmu_415": 0.1188 * (0.19 + 0.0008),
        "b0_dststmu_10413": 0.2727 * 0.1829,
        "b0_dststmu_20413": 0.3069 * 0.4547,
        
        "b+_dststmu_425": 0.1471 * (0.4750 + 0.0014),
        "b+_dststmu_10423": 0.2341 * 0.3821,
        "b+_dststmu_20423": 0.2647 * 0.3821,
    },
    'nbsDstst' : {
        "bs_dststmu_10433": 0.0070 * 0.5,
        "bs_dststmu_435": 0.0070 * 0.5,

    },
    'nDstst_higher' : {
        "b+_dststmu_heavy": 0.05 * 0.08,
        "b0_dststmu_heavy": 0.05 * 0.08,
    },
}

gen_yields = {
    'nDstDX' : {
        "b+_dstdx":1600224718 + 1595754180 ,
        "b0_dstdx": 1256782128 + 1271630284,
    },
    'nbDstst' : {
        "b0_dststmu_415": 1199167539 + 1236752126,
        "b0_dststmu_10413": 1199167539 + 1236752126,
        "b0_dststmu_20413": 1199167539 + 1236752126,
        
        "b+_dststmu_425": 1227691430 + 1175883049,
        "b+_dststmu_10423": 1227691430 + 1175883049,
        "b+_dststmu_20423": 1227691430 + 1175883049,
    },
    'nbsDstst' : {
        "bs_dststmu_10433": 1492068605 + 1453243128,
        "bs_dststmu_435": 1492068605 + 1453243128,

    },
    'nDstst_higher' : {
        "b+_dststmu_heavy": 1338123297 + 1448221903,
        "b0_dststmu_heavy": 1188664181 + 1198834296,
    },
}
bin_labels = {}

if year != 'Run2':
    for i in range(binno):
        print(f'xxxxxxx BIN {i} xxxxxxx')
        for sample in template_label.keys():
            yields_u = []
            errs_u = []
            yields_d = []
            errs_d = []
            with open(fit_setup.template_paths[dim] + f'{year}/' + f'{template_label[sample]}_up_bin{i}{mmiss}_raw.dat', "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    v, e = line.split()
                    yields_u.append(float(v))
                    errs_u.append(float(e))
            with open(fit_setup.template_paths[dim] + f'{year}/' + f'{template_label[sample]}_down_bin{i}{mmiss}_raw.dat', "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    v, e = line.split()
                    yields_d.append(float(v))
                    errs_d.append(float(e))
            h = np.array(yields_u) + np.array(yields_d)
            herr = np.sqrt(np.array(errs_u)**2 + np.array(errs_d)**2)
 
            f = open(fit_setup.template_paths[dim] + f'{year}/' + f'{template_label[sample]}_bin{i}{mmiss}_raw.dat', 'w')
            for j in range(len(h)):

                a = str(h[j]) + ' '
                a += str(herr[j])
                f.write(a+"\n")
            f.close() 
            print(sample, ': ', np.sum(h))



else:
    for i in range(binno):
        print(f'xxxxxxx BIN {i} xxxxxxx')
        for sample in template_label:
            templates = []
            templates_e = []
            for y in ['2016', '2017', '2018']:
                yields = []
                errs = []
                with open(fit_setup.template_paths[dim] + f'/{y}/' + f'{template_label[sample]}_bin{i}{mmiss}_raw.dat', "r") as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        v, e = line.split()
                        yields.append(float(v))
                        errs.append(float(e))
                templates.append(np.array(yields))
                templates_e.append(np.array(errs))
            h = sum(templates)
            print(sample, ' : ', np.sum(h))
            herr = np.sqrt(templates_e[0]**2 + templates_e[1]**2 + templates_e[2]**2 )
            f = open(fit_setup.template_paths[dim] + f'/Run2/' + f'{template_label[sample]}_bin{i}{mmiss}_raw.dat', 'w')
            for j in range(len(h)):

                a = str(h[j]) + ' '
                a += str(herr[j])
                f.write(a+"\n")
            f.close() 
