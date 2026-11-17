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




dim = '1D'
sample = str(sys.argv[1])
year = str(sys.argv[2])
pol = str(sys.argv[3])
template_path = fit_setup.template_path
postfit_bins = 20
postfit_edges = {
    'mu_PE' : np.linspace(0,100000, postfit_bins),
    #'m_corr_reco' : np.linspace(2000,9000, postfit_bins),
    # 'b0_P_Pred' : np.linspace(0, 400000, postfit_bins),
    #'q2GeV_reco_regreg' : np.linspace(0, 11, 5),
}
add_vars = [v for v in postfit_edges.keys()]
json_dir = fit_setup.jsons_dir
bin_vars= fit_setup.angular_vars
bin_file = fit_setup.bin_file
nbins = fit_setup.nbins_4D
ranges = fit_setup.ranges_4D
bins_w = fit_setup.w_bins_4D
bins = json.load(open(json_dir + bin_file))


v0 = bin_vars[0]
v1 = bin_vars[1]
v2 = bin_vars[2]
v3 = bin_vars[3]
print(bins)
print(bin_vars)

WSMu_collapsed_vars = ['ctd_reco_regreg', 'chi_reco_regreg']

def th3d_to_arr(hist):
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    nz = hist.GetNbinsZ()

    values = np.zeros((nx, ny, nz))

    for ix in range(1, nx+1):
        for iy in range(1, ny+1):
            for iz in range(1, nz+1):
                values[ix-1, iy-1, iz-1] = hist.GetBinContent(ix, iy, iz)
    x_edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nx+2)])
    y_edges = np.array([hist.GetYaxis().GetBinLowEdge(i) for i in range(1, ny+2)])
    z_edges = np.array([hist.GetZaxis().GetBinLowEdge(i) for i in range(1, nz+2)])

    return np.nan_to_num(values), (x_edges, y_edges, z_edges)



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



if sample in tuples['MC'].keys():
    data_type = 'MC'

else:
    data_type = 'data'


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
    #"PIDCalibEff",
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
        "b+_dstdx": 'abs(b0_TRUEID) == 521 and abs(dst_TRUEID) == 413',
        "b0_dstdx": 'abs(b0_TRUEID) == 511 and abs(dst_TRUEID) == 413',
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
        "b+_dststmu_heavy": 'abs(b0_TRUEID) == 521 and abs(mu_MC_MOTHER_ID) == 521',
        "b0_dststmu_heavy": 'abs(b0_TRUEID) == 511 and abs(mu_MC_MOTHER_ID) == 511',
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



if sample in tuples['MC'].keys():
    data_type = 'MC'

else:
    data_type = 'data'


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



mc_weight_base = [
    #"PIDCalibEff",
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


r.gInterpreter.Declare("""
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"

double ComputeM2(TVector3 B_0_Dir,
                 TLorentzVector Dst,
                 TLorentzVector D_0,
                 TLorentzVector mu,
                 double pmag) {

    TLorentzVector Vis;

    Vis = Dst + mu;
    double Eb = sqrt(5279.61*5279.61 + pmag*pmag);


    TLorentzVector B_0(
                    pmag * B_0_Dir.Unit().x(),
                    pmag * B_0_Dir.Unit().y(),
                    pmag * B_0_Dir.Unit().z(),
                    Eb
    );
    double m2 = (B_0 - Vis).M2();

    return m2 / 1e6; // GeV^2
}               
                       
double Computew(TVector3 B_0_Dir,
                 TLorentzVector Dst,
                 TLorentzVector D_0,
                 TLorentzVector mu) {

    TLorentzVector B_0, Vis;

    Vis = Dst + mu;

    double scale = 5279.61/Vis.M()*Vis.Z()*B_0_Dir.Mag()/B_0_Dir.Z();
    TVector3 B_vec = scale * B_0_Dir.Unit();

    B_0.SetVectM(B_vec, 5279.61);

    double q2 = (B_0 - Dst).M2();
    double w = (5279.61*5279.61 + Dst.M()*Dst.M() - q2) / (2*5279.61*Dst.M());
    return w;
}
                          
TVector3 RotateVector(TVector3 vec, double angle, TVector3 axis) {
            vec.Rotate(angle, axis);
            return vec;
        }
""")

mmiss_lab = {
    'min' : '_min',
    'max' : '_max',
    'final' : '_final',
    'nom' : '',
    'linear_w_p' : '_linear_p',
    'linear_w_m' : '_linear_m',
    'quad_w_p' : '_quad_p',
    'quad_w_m' : '_quad_m'
}

file_K = r.TFile.Open(f"/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/PID_DDX/effhists-Turbo{year.replace('20','')}-{pol}-K-Brunel_MC15TuneV1_ProbNNk0.2-Brunel_P.Brunel_ETA.nTracks_Brunel.root")
hist_K = file_K.Get("eff_Brunel_MC15TuneV1_ProbNNk>0.2")
file_Pi = r.TFile.Open(f"/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/PID_DDX/effhists-Turbo{year.replace('20','')}-{pol}-Pi-Brunel_MC15TuneV1_ProbNNk0.2-Brunel_P.Brunel_ETA.nTracks_Brunel.root")
hist_Pi = file_Pi.Get("eff_Brunel_MC15TuneV1_ProbNNk>0.2")
file_P = r.TFile.Open(f"/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/PID_DDX/effhists-Turbo{year.replace('20','')}-{pol}-P-Brunel_MC15TuneV1_ProbNNk0.2-Brunel_P.Brunel_ETA.nTracks_Brunel.root")
hist_P = file_P.Get("eff_Brunel_MC15TuneV1_ProbNNk>0.2")
file_Mu = r.TFile.Open(f"/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/PID_DDX/effhists-Turbo{year.replace('20','')}-{pol}-Mu-Brunel_MC15TuneV1_ProbNNk0.2-Brunel_P.Brunel_ETA.nTracks_Brunel.root")
hist_Mu = file_Mu.Get("eff_Brunel_MC15TuneV1_ProbNNk>0.2")
file_e = r.TFile.Open(f"/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/PID_DDX/effhists-Electron{year.replace('20','')}-{pol}-e_B_Jpsi-Brunel_MC15TuneV1_ProbNNk0.2-Brunel_P.Brunel_ETA.nTracks.root")
hist_e = file_e.Get("eff_Brunel_MC15TuneV1_ProbNNk>0.2")
w_K = th3d_to_arr(hist_K)
w_Pi = th3d_to_arr(hist_Pi)
w_P = th3d_to_arr(hist_P)
w_Mu = th3d_to_arr(hist_Mu)
w_e = th3d_to_arr(hist_e)

part_IDs = {

    'e' : 11,
    'Mu' : 13,
    'K' : 321,
    'Pi' : 211,
    'P' : 2212,
}

tables = {
        'e' : w_e,
    'Mu' : w_Mu,
    'K' : w_K,
    'Pi' : w_Pi,
    'P' : w_P,
}

def assign_weights(track_n, df, table):

    weights = np.zeros(len(df))


    for cat, (values, edges) in table.items():
        mask = np.abs(df[f"b0_ISOLATION_TRUEID{track_n}"]) == part_IDs[cat]
        
        if not np.any(mask):
            continue

        x_edges, y_edges, z_edges = edges
        
        ix = np.digitize(df.loc[mask, f"b0_ISOLATION_P{track_n}"], x_edges) - 1
        iy = np.digitize(df.loc[mask, f"b0_ISOLATION_ETA{track_n}"], y_edges) - 1
        iz = np.digitize(df.loc[mask, "nTracks"], z_edges) - 1

        # valid = (
        #     (ix >= 0) & (ix < values.shape[0]) &
        #     (iy >= 0) & (iy < values.shape[1]) &
        #     (iz >= 0) & (iz < values.shape[2])
        # )
        ix = np.clip(ix, 0, len(x_edges) - 2)
        iy = np.clip(iy, 0, len(y_edges) - 2)
        iz = np.clip(iz, 0, len(z_edges) - 2)

        w = np.zeros(np.sum(mask))
        #w[valid] = values[ix[valid], iy[valid], iz[valid]]
        w = values[ix, iy,iz]

        weights[mask] = w

    df[f"PID_weight{track_n}"] = weights
    return df
 


if year == 'Run2':
    for i in range(len(fit_setup.w_edges[5])-1):
        templates = []
        templates_e = []
        for y in ['2016', '2017', '2018']:
            yields = []
            errs = []
            with open(fit_setup.template_path + f'/{y}/' + f'{template_label[sample]}_{pol}_bin{i}_raw.dat', "r") as f:
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
        herr = np.sqrt(templates_e[0]**2 + templates_e[1]**2 + templates_e[2]**2 )
        hn = h / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges['control'][i]))
        herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges['control'][i]))
        f = open(fit_setup.template_path + f'/Run2/' + f'{template_label[sample]}_{pol}_bin{i}.dat', 'w')
        if sample != 'data':
            for j in range(len(hn)):

                a = str(hn[j]) + ' '
                a += str(herrn[j])
                f.write(a+"\n")
            f.close()   
        else:
            for j in range(len(h)):

                a = str(h[j]) + ' '
                a += str(herr[j])
                f.write(a+"\n")
            f.close()  
        f = open(fit_setup.template_path + f'/Run2/' + f'{template_label[sample]}_{pol}_bin{i}_raw.dat', 'w')
        for j in range(len(h)):

            a = str(h[j]) + ' '
            a += str(herr[j])
            f.write(a+"\n")
        f.close() 

else:

    if sample == 'nDstMu' or  sample == 'nDstTau':

        tuple_path = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample]}_{pol}.root'
        tuple_friend = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample]}_{pol}_selection_flags.root'
        tuple_friend2 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample]}_{pol}_addvars.root'
        tuple_friend3 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample]}_{pol}_nspdcorr.root'
        tuple_friend4 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample]}_{pol}_trigemu.root'
        tuple_friend5 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample]}_{pol}_hammer.root'
        tuple_friend6 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample]}_{pol}_gpflow.root'
        tuple_friend7 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample]}_{pol}_wbkin.root'
        tuple_friend8 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample]}_{pol}_angles.root'
        tuple_friend9 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample]}_{pol}_trackeff.root'
        #df_in = r.RDataFrame("DecayTree", str(tuple_path + tuples["MC"][sample]))

        f_main = r.TFile.Open(tuple_path)
        f_friend = r.TFile.Open(tuple_friend)
        f_friend2 = r.TFile.Open(tuple_friend2)
        f_friend3 = r.TFile.Open(tuple_friend3)
        f_friend4 = r.TFile.Open(tuple_friend4)
        f_friend5 = r.TFile.Open(tuple_friend5)
        f_friend6 = r.TFile.Open(tuple_friend6)
        f_friend7 = r.TFile.Open(tuple_friend7)
        f_friend8 = r.TFile.Open(tuple_friend8)
        f_friend9 = r.TFile.Open(tuple_friend9)

        t_main = f_main.Get("DecayTree")
        t_friend = f_friend.Get("DecayTree")
        t_friend2 = f_friend2.Get("DecayTree")
        t_friend3 = f_friend3.Get("DecayTree")
        t_friend4 = f_friend4.Get("DecayTree")
        t_friend5 = f_friend5.Get("DecayTree")
        t_friend6 = f_friend6.Get("DecayTree")
        t_friend7 = f_friend7.Get("DecayTree")
        t_friend8 = f_friend8.Get("DecayTree")
        t_friend9 = f_friend9.Get("DecayTree")
        t_main.AddFriend(t_friend, "sel")
        t_main.AddFriend(t_friend2, "sel2")
        t_main.AddFriend(t_friend3, "sel3")
        t_main.AddFriend(t_friend4, "sel4")
        t_main.AddFriend(t_friend5, "sel5")
        t_main.AddFriend(t_friend6, "sel6")
        t_main.AddFriend(t_friend7, "sel7")
        t_main.AddFriend(t_friend8, "sel8")
        t_main.AddFriend(t_friend9, "sel9")

        df_in = r.RDataFrame(t_main)
        df_sel = df_in.Filter("pass_B0_Kinematic == 1 && pass_K_Kinematic ==1 && pass_Mu_Kinematic ==1 && pass_Pi_Kinematic ==1 && pass_HLT1_Emul ==1  && q2GeV_reco_regreg > 0")#&& pass_PID_Weights ==1
        
        df_sel = df_sel.Define("B_0_Dir", "TVector3(b0_ENDVERTEX_X-b0_OWNPV_X,b0_ENDVERTEX_Y-b0_OWNPV_Y,b0_ENDVERTEX_Z-b0_OWNPV_Z)") \
            .Define("B_0_Dir_TRUE", "TVector3(b0_TRUEP_X, b0_TRUEP_Y, b0_TRUEP_Z)") \
            .Define("B_0_Theta", "abs(B_0_Dir_TRUE.Angle(B_0_Dir))") \
            .Define("B_0_DeltaTheta_Direction", "(B_0_Dir_TRUE.Cross(B_0_Dir)).Unit()") \
            .Define("D0","TLorentzVector(d0_PX,d0_PY,d0_PZ,d0_PE)") \
            .Define("spi", "TLorentzVector(spi_PX,spi_PY,spi_PZ,spi_PE)") \
            .Define("Dst", "D0+spi") \
            .Define("mu", "TLorentzVector(mu_PX,mu_PY,mu_PZ,mu_PE)") \
            .Define("B_0_DeltaTheta_Min","exp(log(B_0_Theta)-0.25) - exp(log(B_0_Theta))") \
            .Define("B_0_DeltaTheta_Max","exp(log(B_0_Theta)+0.25) - exp(log(B_0_Theta))") \
            .Define("B_0_DeltaTheta_Final","exp(log(B_0_Theta)+(0.25*0.35)) - exp(log(B_0_Theta))") \


        df_sel = df_sel.Define(
            "B_0_Dir_Min",
            "RotateVector(B_0_Dir, B_0_DeltaTheta_Min, B_0_DeltaTheta_Direction)"
        )

        df_sel = df_sel.Define(
            "B_0_Dir_Max",
            "RotateVector(B_0_Dir, B_0_DeltaTheta_Max, B_0_DeltaTheta_Direction)"
        )

        df_sel = df_sel.Define(
            "B_0_Dir_Final",
            "RotateVector(B_0_Dir, B_0_DeltaTheta_Final, B_0_DeltaTheta_Direction)"
        )
        df_sel = (df_sel
        .Define("m2_nom", "ComputeM2(B_0_Dir, Dst, D0, mu, B0_reco_approx_P)")
        .Define("m2_min", "ComputeM2(B_0_Dir, Dst, D0, mu, B0_reco_approx_P)")
        .Define("m2_max", "ComputeM2(B_0_Dir, Dst, D0, mu, B0_reco_approx_P)")
        .Define("m2_final", "ComputeM2(B_0_Dir, Dst, D0, mu, B0_reco_approx_P)")
    )
        branch_list = ['w_reco_regreg', 'b0_ISOLATION_BDT', 'b0_ISOLATION_BDT','b0_ISOLATION_BDT2','b0_ISOLATION_BDT3','m2miss_reco_approx', 'nTracks', 'm2_min', 'm2_max', 'm2_nom', 'm2_final', 'B_0_Theta']
        branch_list +=[ 'b0_TRUEID', 'mu_MC_MOTHER_ID', 'dst_TRUEID', 'b0_TrueHadron_D0_ID', 'd0_M', 'dst_M']
        branch_list += mc_weight_base
        branch_list += add_vars
        for coord in ['PX', 'PY', 'PZ', 'TRUEID', 'GhostProb', 'Type']:
            branch_list.append(f'b0_ISOLATION_{coord}')
            branch_list.append(f'b0_ISOLATION_{coord}2')
            branch_list.append(f'b0_ISOLATION_{coord}3')
        if sample in weights_hammer.keys():
            branch_list += [weights_hammer[sample]]
        df = pd.DataFrame(df_sel.AsNumpy(branch_list))
        df['delta_M'] = df['dst_M'] - df['d0_M']
        df = df.query('delta_M >= 140 and delta_M <= 154')
        df_full = len(df)
        print(df_full)

        df = df.query(tm_strings[sample])
        df = df.query('b0_ISOLATION_BDT > 0.3')
        for trck in ['', '2', '3']:
            df[f'b0_ISOLATION_P{trck}'] = np.sqrt(df[f'b0_ISOLATION_PX{trck}']**2 + df[f'b0_ISOLATION_PY{trck}']**2 + df[f'b0_ISOLATION_PZ{trck}']**2)
            df[f'b0_ISOLATION_PT{trck}'] = np.sqrt(df[f'b0_ISOLATION_PX{trck}']**2 + df[f'b0_ISOLATION_PY{trck}']**2)
            df[f'b0_ISOLATION_ETA{trck}'] = df.apply(lambda row: np.arctanh(row[f"b0_ISOLATION_PZ{trck}"] / row[f"b0_ISOLATION_P{trck}"]), axis=1)
        df = assign_weights('', df, tables)
        df = assign_weights('2', df, tables)
        df = assign_weights('3', df, tables)

        df["PID_weight"] = df["PID_weight"] * ((df["b0_ISOLATION_BDT"] > -1.1) & (df["b0_ISOLATION_Type"] == 3))
        df["PID_weight2"] = df["PID_weight2"] * ((df["b0_ISOLATION_BDT2"] > -1.1) & (df["b0_ISOLATION_Type2"] == 3))
        df["PID_weight3"] = df["PID_weight3"] * ((df["b0_ISOLATION_BDT3"] > -1.1) & (df["b0_ISOLATION_Type3"] == 3))
        df["PID_weight_tot"] = df["PID_weight"] + df["PID_weight2"] + df["PID_weight3"] - df["PID_weight"] * df["PID_weight2"] - df["PID_weight"] * df["PID_weight3"] - df["PID_weight2"] * df["PID_weight3"] + df["PID_weight"] * df["PID_weight2"] * df["PID_weight3"]
 
        plt.hist(df["PID_weight_tot"], bins = 100, label = 'PID tot')
        # plt.hist(df['m2_final'], bins = 100, range = [-5,6], histtype='step', label = 'final')
        plt.legend()
        plt.savefig(f'pid_tot_weight_{sample}_{year}_{pol}.pdf')
        print(df["PID_weight_tot"].max())
        print(len(df.query('PID_weight_tot > 100'))/len(df))
        print(df.query('PID_weight_tot > 100')[["b0_ISOLATION_TRUEID","b0_ISOLATION_TRUEID2", 'b0_ISOLATION_TRUEID3', 'b0_ISOLATION_P', 'b0_ISOLATION_ETA', 'nTracks']].head(2))
        df = df.query('PID_weight_tot < 100')
        if sample in weights_hammer.keys():
            df['w_tot'] = df[mc_weight_base + [weights_hammer[sample]] + ['PID_weight_tot']].prod(axis=1)
        else:
            df['w_tot'] = df[mc_weight_base +  ['PID_weight_tot']].prod(axis=1)
        #print(len(df)/df_full)
        df = df.query('w_tot >= 0')
        for mmm in ['m2_min', 'm2_nom', 'm2_max', 'm2miss_reco_approx']:
            print(mmm + ' ', df[mmm].min(), ' ', df[mmm].max())
        for l in ['min', 'nom', 'max', 'final']:
            mmiss = 'm2_' + l
            w_smeared = 'w_' + l
            for i in range(len(fit_setup.w_edges[5])-1):
                
                df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[5][i]} and w_reco_regreg < {fit_setup.w_edges[5][i+1]} and {mmiss} >= {fit_setup.mmiss2_bin_edges["control"][i][0]} and {mmiss} <= {fit_setup.mmiss2_bin_edges["control"][i][-1]}')
                h, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][i], weights = df_w['w_tot'])
                herr = np.sqrt( np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][i], weights=df_w['w_tot']**2, density = False)[0])
                hn, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][i], weights = df_w['w_tot'], density = True)
                herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges['control'][i]))
                f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}_bin{i}{mmiss_lab[l]}.dat', 'w')
                for j in range(len(hn)):

                    a = str(hn[j]) + ' '
                    a += str(herrn[j])
                    f.write(a+"\n")
                f.close()
                f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}_bin{i}{mmiss_lab[l]}_raw.dat', 'w')
                for j in range(len(h)):

                    a = str(h[j]) + ' '
                    a += str(herr[j])
                    f.write(a+"\n")
                f.close()

            df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[5][0]} and w_reco_regreg < {fit_setup.w_edges[5][-1]}')
            h, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges["control"][0], weights = df_w['w_tot'])
            herr = np.sqrt( np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges["control"][0], weights=df_w["w_tot"]**2, density = False)[0])
            hn, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges["control"][0], weights = df_w["w_tot"], density = True)
            herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges["control"][0]))
            f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}{mmiss_lab[l]}_.dat', 'w')
            for j in range(len(hn)):
                a = str(hn[j]) + ' '
                a += str(herrn[j])
                f.write(a+"\n")
            f.close()
            f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}{mmiss_lab[l]}_raw.dat', 'w')
            for j in range(len(h)):
                a = str(h[j]) + ' '
                a += str(herr[j])
                f.write(a+"\n")
            f.close()
        
        print(np.sum(h))

    elif sample in tuples['data'].keys(): 
        if sample in ['data','nWSMu']:
            tuple_path = f'/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/sWeighted_data/{year}/{tuples["data"][sample]}_{pol}_selected_control_withangles_new.root'
            f_main = r.TFile.Open(tuple_path)
            t_main = f_main.Get("DecayTree")
            tuple_friend4 = f'/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/sWeighted_data/{year}/{tuples["data"][sample]}_{pol}_selected_control_sw_withangles_new.root'
            f_friend4 = r.TFile.Open(tuple_friend4)
            t_friend4 = f_friend4.Get("DecayTree")
            t_main.AddFriend(t_friend4, "sel4")
        else:
            tuple_path = fit_setup.tuple_path + f'{year}/{tuples["data"][sample]}_{pol}.root'
            tuple_friend = fit_setup.tuple_path + f'{year}/{tuples["data"][sample]}_{pol}_selection_flags.root'
            tuple_friend2 = fit_setup.tuple_path + f'{year}/{tuples["data"][sample]}_{pol}_addvars.root'
            tuple_friend3 = fit_setup.tuple_path + f'{year}/{tuples["data"][sample]}_{pol}_angles.root'

            #df_in = r.RDataFrame("DecayTree", str(tuple_path + tuples["MC"][sample]))

            f_main = r.TFile.Open(tuple_path)
            f_friend = r.TFile.Open(tuple_friend)
            f_friend2 = r.TFile.Open(tuple_friend2)
            f_friend3 = r.TFile.Open(tuple_friend3)


            t_main = f_main.Get("DecayTree")
            t_friend = f_friend.Get("DecayTree")
            t_friend2 = f_friend2.Get("DecayTree")
            t_friend3 = f_friend3.Get("DecayTree")

            t_main.AddFriend(t_friend, "sel")
            t_main.AddFriend(t_friend2, "sel2")
            t_main.AddFriend(t_friend3, "sel3")

        df_in = r.RDataFrame(t_main)
        if sample in ['data','nWSMu']:
            df_sel = df_in
        else:
            df_sel = df_in.Filter("pass_B0_Kinematic == 1 && pass_K_Kinematic ==1 && pass_Mu_Kinematic ==1 && pass_Pi_Kinematic ==1 && pass_HLT1 ==1 && pass_HLT2 ==1 && pass_nSPD_cut == 1 && pass_L0_Had && q2GeV_reco_regreg > 0 ")#pass_K_PID && pass_Mu_PID && pass_L0_Had && pass_Pi_PID")

        if sample not in ['data','nWSMu']:
            branch_list = ['w_reco_regreg', 'b0_ISOLATION_BDT', 'm2miss_reco_approx', 'b0_ISOLATION_NNk3','b0_ISOLATION_NNk2','b0_ISOLATION_NNk', 'd0_M', 'dst_M']
            branch_list += add_vars
            df = pd.DataFrame(df_sel.AsNumpy(branch_list))
            df['delta_M'] = df['dst_M'] - df['d0_M']
            df = df.query('delta_M >= 140 and delta_M <= 154')
            df = df.query('b0_ISOLATION_BDT > 0.3')
            df['max_prob'] = df[['b0_ISOLATION_NNk', 'b0_ISOLATION_NNk2','b0_ISOLATION_NNk3']].values.max(1)
            df = df.query('max_prob > 0.2')
        else:
            branch_list = ['ns_sw', 'm2miss_reco_approx', 'w_reco_regreg']
            df = pd.DataFrame(df_sel.AsNumpy(branch_list))
        for i in range(len(fit_setup.w_edges[5])-1):
            df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[5][i]} and w_reco_regreg < {fit_setup.w_edges[5][i+1]} and m2miss_reco_approx >= {fit_setup.mmiss2_bin_edges["control"][i][0]} and m2miss_reco_approx <= {fit_setup.mmiss2_bin_edges["control"][i][-1]}')
            if sample in ['data','nWSMu']:
                h, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges['control'][i], weights = df_w['ns_sw'])
                herr = np.sqrt(h)
                hn, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges['control'][i], density = True,  weights = df_w['ns_sw'])
                herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges['control'][i]))
            else:
                h, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges['control'][i])
                herr = np.sqrt(h)
                hn, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges['control'][i], density = True)
                herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges['control'][i]))
            f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}_bin{i}.dat', 'w')
            for j in range(len(h)):
                if sample == 'data':
                    a = str(h[j]) + ' '
                    a += str(herr[j])
                    f.write(a+"\n")
                else:
                    a = str(hn[j]) + ' '
                    a += str(herrn[j])
                    f.write(a+"\n")
            f.close()
            f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}_bin{i}_raw.dat', 'w')
            for j in range(len(h)):
                a = str(h[j]) + ' '
                a += str(herr[j])
                f.write(a+"\n")
            f.close()
            del df_w

        df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[5][0]} and w_reco_regreg < {fit_setup.w_edges[5][-1]}')
        if sample in ['data','nWSMu']:
            h, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges['control'][0], weights = df_w['ns_sw'])
            herr = np.sqrt(h)
            hn, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges['control'][0], density = True, weights = df_w['ns_sw'])
            herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges['control'][0]))
        else:
            h, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges['control'][0])
            herr = np.sqrt(h)
            hn, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges['control'][0], density = True)
            herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges['control'][0]))
        f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}_.dat', 'w')
  
        for j in range(len(h)):
            if sample == 'data':
                a = str(h[j]) + ' '
                a += str(herr[j])
                f.write(a+"\n")
            else:
                a = str(hn[j]) + ' '
                a += str(herrn[j])
                f.write(a+"\n")
        f.close()

        f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}_raw.dat', 'w')
        for j in range(len(h)):
            a = str(h[j]) + ' '
            a += str(herr[j])
            f.write(a+"\n")
        f.close()
        print(np.sum(h))



    else:
        subtemplates = {}
        subtemplates_additional = {}

        for subsample in tm_strings[sample]:
            subtemplates[subsample] = {}
            subtemplates_additional[subsample] = {}
            for v in postfit_edges:
                subtemplates_additional[subsample][v] = {}

            tuple_path = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample][subsample]}_{pol}.root'
            tuple_friend = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample][subsample]}_{pol}_selection_flags.root'
            tuple_friend2 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample][subsample]}_{pol}_addvars.root'
            tuple_friend3 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample][subsample]}_{pol}_nspdcorr.root'
            tuple_friend4 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample][subsample]}_{pol}_trigemu.root'
            tuple_friend5 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample][subsample]}_{pol}_angles.root'
            tuple_friend6 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample][subsample]}_{pol}_gpflow.root'
            tuple_friend7 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample][subsample]}_{pol}_wbkin.root'
            if sample in weights_hammer.keys():
                tuple_friend8 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample][subsample]}_{pol}_hammer.root'
            tuple_friend9 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample][subsample]}_{pol}_trackeff.root'
            #df_in = r.RDataFrame("DecayTree", str(tuple_path + tuples["MC"][sample]))

            f_main = r.TFile.Open(tuple_path)
            f_friend = r.TFile.Open(tuple_friend)
            f_friend2 = r.TFile.Open(tuple_friend2)
            f_friend3 = r.TFile.Open(tuple_friend3)
            f_friend4 = r.TFile.Open(tuple_friend4)
            f_friend5 = r.TFile.Open(tuple_friend5)
            f_friend6 = r.TFile.Open(tuple_friend6)
            f_friend7 = r.TFile.Open(tuple_friend7)
            if sample in weights_hammer.keys():
                f_friend8 = r.TFile.Open(tuple_friend8)
            f_friend9 = r.TFile.Open(tuple_friend9)

            t_main = f_main.Get("DecayTree")
            t_friend = f_friend.Get("DecayTree")
            t_friend2 = f_friend2.Get("DecayTree")
            t_friend3 = f_friend3.Get("DecayTree")
            t_friend4 = f_friend4.Get("DecayTree")
            t_friend5 = f_friend5.Get("DecayTree")
            t_friend6 = f_friend6.Get("DecayTree")
            t_friend7 = f_friend7.Get("DecayTree")
            if sample in weights_hammer.keys():
                t_friend8 = f_friend8.Get("DecayTree")
            t_friend9 = f_friend9.Get("DecayTree")
            t_main.AddFriend(t_friend, "sel")
            t_main.AddFriend(t_friend2, "sel2")
            t_main.AddFriend(t_friend3, "sel3")
            t_main.AddFriend(t_friend4, "sel4")
            t_main.AddFriend(t_friend5, "sel5")
            t_main.AddFriend(t_friend6, "sel6")
            t_main.AddFriend(t_friend7, "sel7")
            if sample in weights_hammer.keys():
                t_main.AddFriend(t_friend8, "sel8")
            t_main.AddFriend(t_friend9, "sel9")

            df_in = r.RDataFrame(t_main)
            df_sel = df_in.Filter("pass_B0_Kinematic == 1 && pass_K_Kinematic ==1 && pass_Mu_Kinematic ==1 && pass_Pi_Kinematic ==1 && pass_HLT1_Emul ==1  && q2GeV_reco_regreg > 0")#&& pass_PID_Weights ==1
            df_sel = df_sel.Define("B_0_Dir", "TVector3(b0_ENDVERTEX_X-b0_OWNPV_X,b0_ENDVERTEX_Y-b0_OWNPV_Y,b0_ENDVERTEX_Z-b0_OWNPV_Z)") \
                .Define("B_0_Dir_TRUE", "TVector3(b0_TRUEP_X, b0_TRUEP_Y, b0_TRUEP_Z)") \
                .Define("B_0_Theta", "abs(B_0_Dir_TRUE.Angle(B_0_Dir))") \
                .Define("B_0_DeltaTheta_Direction", "(B_0_Dir_TRUE.Cross(B_0_Dir)).Unit()") \
                .Define("D0","TLorentzVector(d0_PX,d0_PY,d0_PZ,d0_PE)") \
                .Define("spi", "TLorentzVector(spi_PX,spi_PY,spi_PZ,spi_PE)") \
                .Define("Dst", "D0+spi") \
                .Define("mu", "TLorentzVector(mu_PX,mu_PY,mu_PZ,mu_PE)") \
                .Define("B_0_DeltaTheta_Min","exp(log(B_0_Theta)-0.25) - exp(log(B_0_Theta))") \
                .Define("B_0_DeltaTheta_Max","exp(log(B_0_Theta)+0.25) - exp(log(B_0_Theta))") \
                .Define("B_0_DeltaTheta_Final","exp(log(B_0_Theta)+(0.25*0.35)) - exp(log(B_0_Theta))") \
            
            df_sel = df_sel.Define(
                "B_0_Dir_Min",
                "RotateVector(B_0_Dir, B_0_DeltaTheta_Min, B_0_DeltaTheta_Direction)"
            )

            df_sel = df_sel.Define(
                "B_0_Dir_Max",
                "RotateVector(B_0_Dir, B_0_DeltaTheta_Max, B_0_DeltaTheta_Direction)"
            )
            df_sel = df_sel.Define(
                "B_0_Dir_Final",
                "RotateVector(B_0_Dir, B_0_DeltaTheta_Final, B_0_DeltaTheta_Direction)"
            )

            df_sel = (df_sel
            .Define("m2_nom", "ComputeM2(B_0_Dir, Dst, D0, mu, B0_reco_approx_P)")
            .Define("m2_min", "ComputeM2(B_0_Dir, Dst, D0, mu, B0_reco_approx_P)")
            .Define("m2_max", "ComputeM2(B_0_Dir, Dst, D0, mu, B0_reco_approx_P)")
            .Define("m2_final", "ComputeM2(B_0_Dir, Dst, D0, mu, B0_reco_approx_P)"))


            branch_list = ['w_reco_regreg', 'b0_ISOLATION_BDT', 'b0_ISOLATION_BDT2','b0_ISOLATION_BDT3', 'm2miss_reco_approx', 'nTracks', 'm2_nom', 'm2_min', 'm2_max', 'm2_final']
            branch_list +=[ 'b0_TRUEID', 'mu_MC_MOTHER_ID', 'dst_TRUEID', 'b0_TrueHadron_D0_ID', 'd0_M', 'dst_M']
            branch_list += ['dst_TRUEP_E', 'dst_TRUEP_X','dst_TRUEP_Y','dst_TRUEP_Z', 'mu_MC_MOTHER_TRUEPE', 'mu_MC_MOTHER_TRUEPX', 'mu_MC_MOTHER_TRUEPY', 'mu_MC_MOTHER_TRUEPZ', 'b0_TRUEP_E', 'b0_TRUEP_X','b0_TRUEP_Y','b0_TRUEP_Z', 'k_TRUEP_E', 'k_TRUEP_X','k_TRUEP_Y','k_TRUEP_Z']
            branch_list += mc_weight_base
            branch_list += add_vars

            for coord in ['PX', 'PY', 'PZ', 'TRUEID', 'GhostProb', 'Type']:
                branch_list.append(f'b0_ISOLATION_{coord}')
                branch_list.append(f'b0_ISOLATION_{coord}2')
                branch_list.append(f'b0_ISOLATION_{coord}3')
            if sample in weights_hammer.keys():
                branch_list += [weights_hammer[sample][subsample]]
            df = pd.DataFrame(df_sel.AsNumpy(branch_list))
            df['delta_M'] = df['dst_M'] - df['d0_M']
            df = df.query('delta_M >= 140 and delta_M <= 154')
            df = df.query('b0_ISOLATION_BDT > 0.3')
            for trck in ['', '2', '3']:
                df[f'b0_ISOLATION_P{trck}'] = np.sqrt(df[f'b0_ISOLATION_PX{trck}']**2 + df[f'b0_ISOLATION_PY{trck}']**2 + df[f'b0_ISOLATION_PZ{trck}']**2)
                df[f'b0_ISOLATION_PT{trck}'] = np.sqrt(df[f'b0_ISOLATION_PX{trck}']**2 + df[f'b0_ISOLATION_PY{trck}']**2)
                df[f'b0_ISOLATION_ETA{trck}'] = df.apply(lambda row: np.arctanh(row[f"b0_ISOLATION_PZ{trck}"] / row[f"b0_ISOLATION_P{trck}"]), axis=1)
            df_full = len(df)
            print(df_full)

            df = df.query(tm_strings[sample][subsample])
            df = assign_weights('', df, tables)
            df = assign_weights('2', df, tables)
            df = assign_weights('3', df, tables)

            df["PID_weight"] = df["PID_weight"] * ((df["b0_ISOLATION_BDT"] > -1.1) & (df["b0_ISOLATION_Type"] == 3))
            df["PID_weight2"] = df["PID_weight2"] * ((df["b0_ISOLATION_BDT2"] > -1.1) & (df["b0_ISOLATION_Type2"] == 3))
            df["PID_weight3"] = df["PID_weight3"] * ((df["b0_ISOLATION_BDT3"] > -1.1) & (df["b0_ISOLATION_Type3"] == 3))
            df["PID_weight_tot"] = df["PID_weight"] + df["PID_weight2"] + df["PID_weight3"] - df["PID_weight"] * df["PID_weight2"] - df["PID_weight"] * df["PID_weight3"] - df["PID_weight2"] * df["PID_weight3"] + df["PID_weight"] * df["PID_weight2"] * df["PID_weight3"]
            
            df['b0_invm'] = np.sqrt(df['b0_TRUEP_E']**2 - (df['b0_TRUEP_X']**2 + df['b0_TRUEP_Y']**2 + df['b0_TRUEP_Z']**2))
            df['dst_invm'] = np.sqrt(df['dst_TRUEP_E']**2 - (df['dst_TRUEP_X']**2 + df['dst_TRUEP_Y']**2 + df['dst_TRUEP_Z']**2))
            df['k_invm'] = np.sqrt(df['k_TRUEP_E']**2 - (df['k_TRUEP_X']**2 + df['k_TRUEP_Y']**2 + df['k_TRUEP_Z']**2))
            df['d_invm'] = np.sqrt(df['mu_MC_MOTHER_TRUEPE']**2 - (df['mu_MC_MOTHER_TRUEPX']**2 + df['mu_MC_MOTHER_TRUEPY']**2 + df['mu_MC_MOTHER_TRUEPZ']**2))
            df['dd_invm'] = np.sqrt((df['dst_TRUEP_E'] + df['mu_MC_MOTHER_TRUEPE'])**2 - ((df['dst_TRUEP_X'] + df['mu_MC_MOTHER_TRUEPX'])**2 + (df['dst_TRUEP_Y'] + df['mu_MC_MOTHER_TRUEPY'])**2 + (df['dst_TRUEP_Z'] + df['mu_MC_MOTHER_TRUEPZ'])**2))


            df['linear_w_p'] = np.where( df['dd_invm'] < 4990, 1 + 2 * (np.sqrt((df['dd_invm']**2 - (df['dst_invm'] + df['d_invm'])**2)/((df['b0_invm'] - df['k_invm'])**2 - (df['dst_invm'] + df['d_invm'])**2)) - 0.5 ) ,1) 
            df['linear_w_m'] = np.where( df['dd_invm'] < 4990,1 - 2 * (np.sqrt((df['dd_invm']**2 - (df['dst_invm'] + df['d_invm'])**2)/((df['b0_invm'] - df['k_invm'])**2 - (df['dst_invm'] + df['d_invm'])**2)) - 0.5 ) ,1)        
            df['quad_w_p'] =  np.where( df['dd_invm'] < 4990,8 * (np.sqrt((df['dd_invm']**2 - (df['dst_invm'] + df['d_invm'])**2)/((df['b0_invm'] - df['k_invm'])**2 - (df['dst_invm'] + df['d_invm'])**2)) - 0.5 )**2 ,1)  
            df['quad_w_m'] = np.where( df['dd_invm'] < 4990,2 - 8 * (np.sqrt((df['dd_invm']**2 - (df['dst_invm'] + df['d_invm'])**2)/((df['b0_invm'] - df['k_invm'])**2 - (df['dst_invm'] + df['d_invm'])**2)) - 0.5 )**2  ,1)       
            
            plt.hist(df['linear_w_m'], bins = 100, histtype='step')
            plt.savefig('linear_m.pdf') 
            plt.clf()
            plt.hist(df['linear_w_p'], bins = 100, histtype='step')
            plt.savefig('linear_p.pdf') 
            plt.clf()
            plt.hist(df['quad_w_m'], bins = 100, histtype='step')
            plt.savefig('quad_m.pdf') 
            plt.clf()
            plt.hist(df['quad_w_p'], bins = 100, histtype='step')
            plt.savefig('quad_p.pdf') 
            plt.clf()
            plt.hist(df['PID_weight_tot'], bins = 100, histtype='step')
            plt.savefig(f'PID_weight_{subsample}_{year}.pdf') 
            plt.clf()
            # print(df.head(10))
            if sample in weights_hammer.keys():
                df['w_tot'] = df[mc_weight_base + [weights_hammer[sample][subsample]] + ['PID_weight_tot']].prod(axis=1)
            else:
                df['w_tot'] = df[mc_weight_base  + ['PID_weight_tot']].prod(axis=1)
            print(subsample, ' pct of neg wgts: ', len(df.query('w_tot < 0'))/len(df))
            #print(len(df)/df_full)
            df = df.query('w_tot >= 0')
            for l in ['min', 'nom', 'max', 'final']:
                mmiss = 'm2_' + l
                w_smeared = 'w_' + l
                subtemplates[subsample][l] = {}
                for i in range(len(fit_setup.w_edges[5])-1):
                    df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[5][i]} and w_reco_regreg < {fit_setup.w_edges[5][i+1]}  and {mmiss} >= {fit_setup.mmiss2_bin_edges["control"][i][0]} and {mmiss} <= {fit_setup.mmiss2_bin_edges["control"][i][-1]}')
                    h, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][i], weights = df_w['w_tot'])
                    hn, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][i], weights = df_w['w_tot'], density = True)
                    herr = np.sqrt( np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][i], weights=df_w['w_tot']**2, density = False)[0])
                    herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges['control'][i]))
                    sel_evts = np.sum(h)
                    subtemplates[subsample][l][i] = {'nev': np.sum(h),'template': hn, 'template_err': herrn, 'contrib': rel_BFs[sample][subsample][0] *(sel_evts) / (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample]), 'contrib_error': np.sqrt((rel_BFs[sample][subsample][0]**2  * sel_evts**2) * ((rel_BFs[sample][subsample][1]/ rel_BFs[sample][subsample][0])**2 + (np.sqrt(sel_evts)/sel_evts)**2) )/ (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample])}

                    del df_w
                df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[5][0]} and w_reco_regreg < {fit_setup.w_edges[5][-1]}  and {mmiss} >= {fit_setup.mmiss2_bin_edges["control"][0][0]} and {mmiss} <= {fit_setup.mmiss2_bin_edges["control"][0][-1]}')
                h, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][0], weights = df_w['w_tot'])
                hn, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][0], weights = df_w['w_tot'], density = True)
                herr = np.sqrt( np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][0], weights=df_w['w_tot']**2, density = False)[0])
                herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges['control'][0]))
                sel_evts = np.sum(h)
                subtemplates[subsample][l]['full'] = {'nev': np.sum(h),'template': hn, 'template_err': herrn, 'contrib': rel_BFs[sample][subsample][0] *(sel_evts) / (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample]), 'contrib_error': np.sqrt((rel_BFs[sample][subsample][0]**2  * sel_evts**2) * ((rel_BFs[sample][subsample][1]/ rel_BFs[sample][subsample][0])**2 + (np.sqrt(sel_evts)/sel_evts)**2) )/ (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample])}

                del df_w
            for l in ['linear_w_p', 'linear_w_m', 'quad_w_p', 'quad_w_m']:
                subtemplates[subsample][l] = {}
                mmiss = 'm2_nom'
                for i in range(len(fit_setup.w_edges[5])-1):
                    df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[5][i]} and w_reco_regreg < {fit_setup.w_edges[5][i+1]}  and {mmiss} >= {fit_setup.mmiss2_bin_edges["control"][i][0]} and {mmiss} <= {fit_setup.mmiss2_bin_edges["control"][i][-1]}')
                    h, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][i], weights = df_w['w_tot']*df_w[l])
                    hn, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][i], weights = df_w['w_tot']*df_w[l], density = True)
                    herr = np.sqrt( np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][i], weights=(df_w['w_tot']*df_w[l])**2, density = False)[0])
                    herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges['control'][i]))
                    sel_evts = np.sum(h)
                    subtemplates[subsample][l][i] = {'nev': np.sum(h),'template': hn, 'template_err': herrn, 'contrib': rel_BFs[sample][subsample][0] *(sel_evts) / (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample]), 'contrib_error': np.sqrt((rel_BFs[sample][subsample][0]**2  * sel_evts**2) * ((rel_BFs[sample][subsample][1]/ rel_BFs[sample][subsample][0])**2 + (np.sqrt(sel_evts)/sel_evts)**2) )/ (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample])}

                    del df_w
                df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[5][0]} and w_reco_regreg < {fit_setup.w_edges[5][-1]}  and {mmiss} >= {fit_setup.mmiss2_bin_edges["control"][0][0]} and {mmiss} <= {fit_setup.mmiss2_bin_edges["control"][0][-1]}')
                h, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][0], weights = df_w['w_tot']*df_w[l])
                hn, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][0], weights = df_w['w_tot']*df_w[l], density = True)
                herr = np.sqrt( np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges['control'][0], weights=(df_w['w_tot']*df_w[l])**2, density = False)[0])
                herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges['control'][0]))
                sel_evts = np.sum(h)
                subtemplates[subsample][l]['full'] = {'nev': np.sum(h),'template': hn, 'template_err': herrn, 'contrib': rel_BFs[sample][subsample][0] *(sel_evts) / (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample]), 'contrib_error': np.sqrt((rel_BFs[sample][subsample][0]**2  * sel_evts**2) * ((rel_BFs[sample][subsample][1]/ rel_BFs[sample][subsample][0])**2 + (np.sqrt(sel_evts)/sel_evts)**2) )/ (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample])}

                del df_w
        for l in ['min', 'nom', 'max', 'final','linear_w_p', 'linear_w_m', 'quad_w_p', 'quad_w_m']:
            for i in range(len(fit_setup.w_edges[5])-1):
                full_temp = np.zeros(len(fit_setup.mmiss2_bin_edges['control'][i])-1)
                full_err = np.zeros(len(fit_setup.mmiss2_bin_edges['control'][i])-1)
                num_err = np.zeros(len(fit_setup.mmiss2_bin_edges['control'][i])-1)
                den_err = np.zeros(len(fit_setup.mmiss2_bin_edges['control'][i])-1)
                full_contribs = [subtemplates[k][l][i]['contrib']  for k in subtemplates.keys()]
                full_contrib = np.sum(full_contribs)

                nevs = [subtemplates[k][l][i]['nev']  for k in subtemplates.keys()]
                nev = np.sum(nevs)
                for subsample in subtemplates.keys():
                    full_temp += subtemplates[subsample][l][i]['contrib'] * subtemplates[subsample][l][i]['template']
                    num_err_ = (subtemplates[subsample][l][i]['contrib'] * subtemplates[subsample][l][i]['template'])**2 * ( (subtemplates[subsample][l][i]['contrib_error']/subtemplates[subsample][l][i]['contrib'])**2 + (subtemplates[subsample][l][i]['template_err']/subtemplates[subsample][l][i]['template'])**2)
                    num_err_ = np.nan_to_num(num_err_)
                    num_err += num_err_
                    den_err += subtemplates[subsample][l][i]['contrib_error']**2

                full_err = (full_temp/full_contrib)**2 * ( (num_err/ full_temp**2)  + (den_err/ full_contrib**2))
                full_err = np.nan_to_num(full_err)
                full_err = np.sqrt(full_err)
                #full_err = np.sqrt(full_err)
                full_temp /= full_contrib

                f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}_bin{i}{mmiss_lab[l]}.dat', 'w')
                for j in range(len(full_temp)):
                    a = str(full_temp[j]) + ' '
                    a += str(full_err[j])
                    f.write(a+"\n")
                f.close()
                full_temp *= nev * np.diff(fit_setup.mmiss2_bin_edges['control'][i])
                full_err *= nev * np.diff(fit_setup.mmiss2_bin_edges['control'][i])
                f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}_bin{i}{mmiss_lab[l]}_raw.dat', 'w')
                for j in range(len(full_temp)):
                    a = str(full_temp[j]) + ' '
                    a += str(full_err[j])
                    f.write(a+"\n")
                f.close()

            full_temp = np.zeros(len(fit_setup.mmiss2_bin_edges['control'][0])-1)
            full_err = np.zeros(len(fit_setup.mmiss2_bin_edges['control'][0])-1)
            num_err = np.zeros(len(fit_setup.mmiss2_bin_edges['control'][0])-1)
            den_err = np.zeros(len(fit_setup.mmiss2_bin_edges['control'][0])-1)
            full_contribs = [subtemplates[k][l]['full']['contrib']  for k in subtemplates.keys()]
            full_contrib = np.sum(full_contribs)

            nevs = [subtemplates[k][l]['full']['nev']  for k in subtemplates.keys()]
            nev = np.sum(nevs)
            for subsample in subtemplates.keys():
                full_temp += subtemplates[subsample][l]['full']['contrib'] * subtemplates[subsample][l]['full']['template']
                num_err_ = (subtemplates[subsample][l]['full']['contrib'] * subtemplates[subsample][l]['full']['template'])**2 * ( (subtemplates[subsample][l]['full']['contrib_error']/subtemplates[subsample][l]['full']['contrib'])**2 + (subtemplates[subsample][l]['full']['template_err']/subtemplates[subsample][l]['full']['template'])**2)
                num_err_ = np.nan_to_num(num_err_)
                num_err += num_err_
                den_err += subtemplates[subsample][l]['full']['contrib_error']**2

            full_err = (full_temp/full_contrib)**2 * ( (num_err/ full_temp**2)  + (den_err/ full_contrib**2))
            full_err = np.nan_to_num(full_err)
            full_err = np.sqrt(full_err)
            #full_err = np.sqrt(full_err)
            full_temp /= full_contrib

            f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}{mmiss_lab[l]}.dat', 'w')
            for j in range(len(full_temp)):
                a = str(full_temp[j]) + ' '
                a += str(full_err[j])
                f.write(a+"\n")
            f.close()
            full_temp *= nev * np.diff(fit_setup.mmiss2_bin_edges['control'][0])
            full_err *= nev * np.diff(fit_setup.mmiss2_bin_edges['control'][0])
            f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}{mmiss_lab[l]}_raw.dat', 'w')
            for j in range(len(full_temp)):
                a = str(full_temp[j]) + ' '
                a += str(full_err[j])
                f.write(a+"\n")
            f.close()


