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




dim = '5D'
n_w_bins = 5

year = str(sys.argv[1])
blind  = True
bfac = 1
if blind:
    bfac = np.random.rand()

template_path = fit_setup.template_paths[dim]
postfit_bins = 20
sample = 'nDstMu'
edges = dim + f'_{n_w_bins}'

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

smear = {
    '2016': 0.2732,
    '2017': 0.3336,
    '2018': 0.2782,
}


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














r.gInterpreter.Declare("""
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <cmath>

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
                       
double Computew(    double pmag,
                       TVector3 B_0_Dir,
                 TLorentzVector Dst) {

    double Eb = sqrt(5279.61 * 5279.61 + pmag*pmag);


    TLorentzVector B_0(
                    pmag * B_0_Dir.Unit().x(),
                    pmag * B_0_Dir.Unit().y(),
                    pmag * B_0_Dir.Unit().z(),
                    Eb
    );

    double q2 = (B_0 - Dst).M2();
    double w = (5279.61*5279.61 + Dst.M()*Dst.M() - q2) / (2*5279.61*Dst.M());
    return w;
}
                          
TVector3 RotateVector(TVector3 vec, double angle, TVector3 axis) {
            vec.Rotate(angle, axis);
            return vec;
        }



                       
double ComputeCosThetaL(
    double pmag,
    TVector3 B_0_Dir,
    TLorentzVector Dst,
    TLorentzVector mu)
{
    TLorentzVector Vis = Dst + mu;

    double Eb = sqrt(5279.61 * 5279.61 + pmag*pmag);


    TLorentzVector B_0(
                    pmag * B_0_Dir.Unit().x(),
                    pmag * B_0_Dir.Unit().y(),
                    pmag * B_0_Dir.Unit().z(),
                    Eb
    );


    TVector3 boostB = -B_0.BoostVector();
       
    TLorentzVector Dst_B = Dst;
    TLorentzVector mu_B  = mu;

    Dst_B.Boost(boostB);
    mu_B.Boost(boostB);

    TLorentzVector W_B(
        -Dst_B.Px(),
        -Dst_B.Py(),
        -Dst_B.Pz(),
        B_0.M() - Dst_B.E()
    );

    if (W_B.M2() <= 0.)
        return -999.;
    TVector3 boostW = -W_B.BoostVector();  
    TLorentzVector mu_W = mu_B;
    mu_W.Boost(boostW);

    return mu_W.Vect().Unit().Dot(
        Dst_B.Vect().Unit()
    );
}
                       
double ComputeCosThetaD(
    double pmag,
    TVector3 B_0_Dir,
    TLorentzVector Dst,
    TLorentzVector D)
{
    double Eb = sqrt(5279.61 * 5279.61 + pmag*pmag);


    TLorentzVector B_0(
                    pmag * B_0_Dir.Unit().x(),
                    pmag * B_0_Dir.Unit().y(),
                    pmag * B_0_Dir.Unit().z(),
                    Eb
    );


    TVector3 boostB = -B_0.BoostVector();
       
    TLorentzVector Dst_B = Dst;
    TLorentzVector D_B  = D;

    Dst_B.Boost(boostB);
    D_B.Boost(boostB);

    TVector3 boostDst = -Dst_B.BoostVector();
    TLorentzVector D_Dst  = D_B;
    D_Dst.Boost(boostDst);
                       

    return D_Dst.Vect().Unit().Dot(
        Dst_B.Vect().Unit()
    );
}
                       
double ComputeCosChi(
    double pmag,
    TVector3 B_0_Dir,
    TLorentzVector Dst,
    TLorentzVector D,
    TLorentzVector Mu,
    TLorentzVector Pi)
{
    double Eb = sqrt(5279.61 * 5279.61 + pmag*pmag);


    TLorentzVector B_0(
                    pmag * B_0_Dir.Unit().x(),
                    pmag * B_0_Dir.Unit().y(),
                    pmag * B_0_Dir.Unit().z(),
                    Eb
    );


    TVector3 boostB = -B_0.BoostVector();
       
    TLorentzVector Dst_B = Dst;
    TLorentzVector D_B  = D;
    TLorentzVector Mu_B  = Mu;
    TLorentzVector Pi_B  = Pi;

    Dst_B.Boost(boostB);
    D_B.Boost(boostB);
    Mu_B.Boost(boostB);
    Pi_B.Boost(boostB);
                       
    TVector3 Dst_B_p3 = Dst_B.Vect().Unit();
    TVector3 D_B_p3 = D_B.Vect().Unit();
    TVector3 Mu_B_p3 = Mu_B.Vect().Unit();
    TVector3 Pi_B_p3 = Pi_B.Vect().Unit();
                       
    TVector3 cross_Mu_Dst = (Mu_B_p3.Cross(-Dst_B_p3)).Unit();
    TVector3 cross_D0_Pi = (D_B_p3.Cross(Pi_B_p3)).Unit();
                       

    return cross_Mu_Dst.Dot(cross_D0_Pi);
    
}
                       
double ComputeSinChi(
    double pmag,
    TVector3 B_0_Dir,
    TLorentzVector Dst,
    TLorentzVector D,
    TLorentzVector Mu,
    TLorentzVector Pi,
    int BID)
{
    double Eb = sqrt(5279.61 * 5279.61 + pmag*pmag);


    TLorentzVector B_0(
                    pmag * B_0_Dir.Unit().x(),
                    pmag * B_0_Dir.Unit().y(),
                    pmag * B_0_Dir.Unit().z(),
                    Eb
    );


    TVector3 boostB = -B_0.BoostVector();
       
    TLorentzVector Dst_B = Dst;
    TLorentzVector D_B  = D;
    TLorentzVector Mu_B  = Mu;
    TLorentzVector Pi_B  = Pi;

    Dst_B.Boost(boostB);
    D_B.Boost(boostB);
    Mu_B.Boost(boostB);
    Pi_B.Boost(boostB);
                       
    TVector3 Dst_B_p3 = Dst_B.Vect().Unit();
    TVector3 D_B_p3 = D_B.Vect().Unit();
    TVector3 Mu_B_p3 = Mu_B.Vect().Unit();
    TVector3 Pi_B_p3 = Pi_B.Vect().Unit();
                       
    TVector3 cross_Mu_Dst = (Mu_B_p3.Cross(-Dst_B_p3)).Unit();
    TVector3 cross_D0_Pi = (D_B_p3.Cross(Pi_B_p3)).Unit();
    TVector3 cross_vecs = (cross_Mu_Dst.Cross(cross_D0_Pi));

    double sinchi = cross_vecs.Dot(Dst_B_p3);              

    return BID >= 0 ? sinchi : -sinchi;
    
}
                       
double ComputeChi(
    double pmag,
    TVector3 B_0_Dir,
    TLorentzVector Dst,
    TLorentzVector D,
    TLorentzVector Mu,
    TLorentzVector Pi,
    int BID)
{
    double coschi = ComputeCosChi(pmag, B_0_Dir, Dst, D, Mu, Pi);
    double sinchi = ComputeSinChi(pmag, B_0_Dir, Dst, D, Mu, Pi, BID);
    double chi = TMath::ACos(coschi);
    return sinchi >= 0 ? chi : -chi;
    
}
""")


mmiss_lab = {
    'min' : '_min',
    'max' : '_max',
    'final' : '_final',
    'nom' : ''
}

bin_labels = {}
binno = 0
for bv0 in range(len(bins[v0])-1):
    for bv1 in range(len(bins[v1])-1):
        for bv2 in range(len(bins[v2])-1):
            for bv3 in range(len(bins[v3])-1):
                bin_labels[binno] = f'{bv0}_{bv1}'
                binno+=1


tuple_path_gen = f'/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/{year}_gen_level_tuple.root'

f_main_gen  = r.TFile.Open(tuple_path_gen )

t_main_gen  = f_main_gen .Get("MCDecayTreeTuple/MCDecayTree")

df_sel_gen  = r.RDataFrame(t_main_gen )

df_sel_gen  = (df_sel_gen 
    .Alias("D0_TRUEP_X", "D~0_TRUEP_X")
    .Alias("D0_TRUEP_Y", "D~0_TRUEP_Y")
    .Alias("D0_TRUEP_Z", "D~0_TRUEP_Z")
    .Alias("D0_TRUEP_E", "D~0_TRUEP_E")
)
df_sel_gen  =  df_sel_gen .Define("B_0_Dir", "TVector3(B0_TRUEP_X, B0_TRUEP_Y, B0_TRUEP_Z)") \
    .Define("B0","TLorentzVector(B0_TRUEP_X,B0_TRUEP_Y,B0_TRUEP_Z,B0_TRUEP_E)") \
    .Define("B0_P","B0.P()") \
    .Define("D0","TLorentzVector(D0_TRUEP_X,D0_TRUEP_Y,D0_TRUEP_Z,D0_TRUEP_E)") \
    .Define("spi", "TLorentzVector(piminus0_TRUEP_X,piminus0_TRUEP_Y,piminus0_TRUEP_Z,piminus0_TRUEP_E)") \
    .Define("Dst", "TLorentzVector(Dst_2010_minus_TRUEP_X,Dst_2010_minus_TRUEP_Y,Dst_2010_minus_TRUEP_Z,Dst_2010_minus_TRUEP_E)") \
    .Define("Dst_sum", "D0+spi") \
    .Define("Dst_sum_m", "Dst_sum.M()") \
    .Define("Dst_m", "Dst.M()") \
    .Define("mu", "TLorentzVector(muplus_TRUEP_X,muplus_TRUEP_Y,muplus_TRUEP_Z,muplus_TRUEP_E)") \



df_sel_gen  = df_sel_gen .Define("ctl_true","ComputeCosThetaL(B0_P,B_0_Dir, Dst, mu)")

df_sel_gen  = df_sel_gen .Define("ctd_true","ComputeCosThetaD(B0_P,B_0_Dir, Dst, D0)")

df_sel_gen  = df_sel_gen .Define("chi_true","ComputeChi(B0_P,B_0_Dir, Dst, D0, mu, spi, Dst_2010_minus_MC_MOTHER_ID)")

df_sel_gen  = df_sel_gen .Define("w_true","Computew(B0_P,B_0_Dir, Dst)")


branch_list = ['w_true', 'ctl_true', 'chi_true', 'ctd_true']#, 'w_min', 'w_max', 'w_nom']



df_gen  = pd.DataFrame(df_sel_gen .AsNumpy(branch_list))
print(v0,v1,v2,v3)
hist, edges = np.histogramdd(
    df_gen[["w_true", "ctl_true", "ctd_true", "chi_true"]].to_numpy(),
    bins=(bins[v0], bins[v1], bins[v2], bins[v3])
)

hist1d = hist.ravel()

yields_gen  = []
if dim == '1D':
    for i in range(len(fit_setup.w_edges[n_w_bins])-1):

        df_w_gen = df_gen .query(f'w_true >= {fit_setup.w_edges[n_w_bins][i]} and w_true < {fit_setup.w_edges[n_w_bins][i+1]}')
        yields_gen .append(len(df_w_gen ))
        del df_w_gen 

elif dim == '5D':

    bin_no = 0

    for iv0, (bv0l, bv0h) in enumerate(zip(bins[v0], bins[v0][1:])):
        for iv1, (bv1l, bv1h )in enumerate(zip(bins[v1], bins[v1][1:])):
            for bv2l, bv2h in zip(bins[v2], bins[v2][1:]):
                for bv3l, bv3h in zip(bins[v3], bins[v3][1:]):
                    v0_ = v0.replace('reco_regreg','true')
                    v1_ = v1.replace('reco_regreg','true')
                    v2_ = v2.replace('reco_regreg','true')
                    v3_ = v3.replace('reco_regreg','true')
                    bin_sel = f'{v0_} >= {bv0l} and {v0_} < {bv0h} '
                    bin_sel += f'and {v1_} >= {bv1l} and {v1_} < {bv1h} '
                    bin_sel += f'and {v2_} >= {bv2l} and {v2_} < {bv2h} '
                    bin_sel += f'and {v3_} >= {bv3l} and {v3_} < {bv3h}'

                    #mmiss_bins = bins['m2miss_reco_approx'][f'{iv0}_{iv1}']
                    #mmiss_bins = fit_setup.mmiss2_bin_edges['1D_5'][iv0]
                    df_w_gen  = df_gen .query(bin_sel)


                    yields_gen .append(len(df_w_gen ))
                    del df_w_gen 
                    bin_no += 1


errs_gen = np.sqrt(np.array(yields_gen))


yields_up_down = {}
errs_up_down = {}
for pol in ['up', 'down']:
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
    tuple_friend10 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample]}_{pol}_pidcalib.root'
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
    f_friend10 = r.TFile.Open(tuple_friend10)

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
    t_friend10 = f_friend10.Get("DecayTree")

    t_main.AddFriend(t_friend, "sel")
    t_main.AddFriend(t_friend2, "sel2")
    t_main.AddFriend(t_friend3, "sel3")
    t_main.AddFriend(t_friend4, "sel4")
    t_main.AddFriend(t_friend5, "sel5")
    t_main.AddFriend(t_friend6, "sel6")
    t_main.AddFriend(t_friend7, "sel7")
    t_main.AddFriend(t_friend8, "sel8")
    t_main.AddFriend(t_friend9, "sel9")
    t_main.AddFriend(t_friend10, "sel10")

    df_in = r.RDataFrame(t_main)
    #df_sel = df_in.Filter("pass_B0_Kinematic == 1 && pass_D0_Kinematic ==1 && pass_Dst_Kinematic ==1 && pass_K_Kinematic ==1 && pass_Mu_Kinematic ==1 && pass_Pi_Kinematic ==1 && pass_HLT1_Emul ==1  && q2GeV_reco_regreg > 0 && pass_PID_Weights ==1 && pass_B0_Isolation == 1")
    df_sel = df_in.Filter("pass_B0_Kinematic == 1  && pass_K_Kinematic ==1 && pass_Mu_Kinematic ==1 && pass_Pi_Kinematic ==1 && pass_HLT1_Emul ==1  && q2GeV_reco_regreg > 0 && pass_PID_Weights ==1 && pass_B0_Isolation == 1 && m2miss_reco_approx >= -5 &&  m2miss_reco_approx <= 6 && w_reco_regreg >= 1 && w_reco_regreg <= 1.51")



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
        .Define("B_0_DeltaTheta_Final",f"exp(log(B_0_Theta)+(0.25*{smear[year]})) - exp(log(B_0_Theta))") \
        .Define("B_0_TRUE_P4", "TLorentzVector(b0_TRUEP_X, b0_TRUEP_Y, b0_TRUEP_Z, b0_TRUEP_E)") \
        .Define("B_0_TRUE_P", "B_0_TRUE_P4.P()") \
        .Define("D0_TRUE","TLorentzVector(d0_TRUEP_X,d0_TRUEP_Y,d0_TRUEP_Z,d0_TRUEP_E)") \
        .Define("spi_TRUE", "TLorentzVector(spi_TRUEP_X,spi_TRUEP_Y,spi_TRUEP_Z,spi_TRUEP_E)") \
        .Define("Dst_TRUE", "D0_TRUE+spi_TRUE") \
        .Define("mu_TRUE", "TLorentzVector(mu_TRUEP_X,mu_TRUEP_Y,mu_TRUEP_Z,mu_TRUEP_E)") \


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
    .Define("m2_regreg", "ComputeM2(B_0_Dir, Dst, D0, mu, B0_reco_regreg_P)")
    .Define("m2_min", "ComputeM2(B_0_Dir_Min, Dst, D0, mu, B0_reco_approx_P)")
    .Define("m2_max", "ComputeM2(B_0_Dir_Max, Dst, D0, mu, B0_reco_approx_P)")
    .Define("m2_final", "ComputeM2(B_0_Dir_Final, Dst, D0, mu, B0_reco_approx_P)")
    )
    df_sel = df_sel.Define("ctl_reco_nom","ComputeCosThetaL(B0_reco_regreg_P,B_0_Dir, Dst, mu)")
    df_sel = df_sel.Define("ctl_reco_min","ComputeCosThetaL(B0_reco_regreg_P,B_0_Dir_Min, Dst, mu)")
    df_sel = df_sel.Define("ctl_reco_max","ComputeCosThetaL(B0_reco_regreg_P,B_0_Dir_Max, Dst, mu)")
    df_sel = df_sel.Define("ctl_reco_final","ComputeCosThetaL(B0_reco_regreg_P,B_0_Dir_Final, Dst, mu)")
    df_sel = df_sel.Define("ctl_true","ComputeCosThetaL(B_0_TRUE_P,B_0_Dir_TRUE, Dst_TRUE, mu_TRUE)")

    df_sel = df_sel.Define("ctd_reco_nom","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir, Dst, D0)")
    df_sel = df_sel.Define("ctd_reco_min","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir_Min, Dst, D0)")
    df_sel = df_sel.Define("ctd_reco_max","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir_Max, Dst, D0)")
    df_sel = df_sel.Define("ctd_reco_final","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir_Final, Dst, D0)")
    df_sel = df_sel.Define("ctd_true","ComputeCosThetaD(B_0_TRUE_P,B_0_Dir_TRUE, Dst_TRUE, D0_TRUE)")

    df_sel = df_sel.Define("chi_reco_nom","ComputeChi(B0_reco_regreg_P,B_0_Dir, Dst, D0, mu, spi, b0_ID)")
    df_sel = df_sel.Define("chi_reco_min","ComputeChi(B0_reco_regreg_P,B_0_Dir_Min, Dst, D0, mu, spi,b0_ID)")
    df_sel = df_sel.Define("chi_reco_max","ComputeChi(B0_reco_regreg_P,B_0_Dir_Max, Dst,mu, spi, D0, b0_ID)")
    df_sel = df_sel.Define("chi_reco_final","ComputeChi(B0_reco_regreg_P,B_0_Dir_Final, Dst, D0,mu, spi, b0_ID)")
    df_sel = df_sel.Define("chi_true","ComputeChi(B_0_TRUE_P,B_0_Dir_TRUE, Dst_TRUE, D0_TRUE,mu_TRUE, spi_TRUE, b0_TRUEID)")


    df_sel = df_sel.Define("w_reco_nom","Computew(B0_reco_regreg_P,B_0_Dir, Dst)")
    df_sel = df_sel.Define("w_reco_min","Computew(B0_reco_regreg_P,B_0_Dir_Min, Dst)")
    df_sel = df_sel.Define("w_reco_max","Computew(B0_reco_regreg_P,B_0_Dir_Max, Dst)")
    df_sel = df_sel.Define("w_reco_final","Computew(B0_reco_regreg_P,B_0_Dir_Final, Dst)")
    df_sel = df_sel.Define("w_true","Computew(B_0_TRUE_P,B_0_Dir_TRUE, Dst_TRUE)")

    branch_list = ['w_reco_regreg', 'm2miss_reco_approx', 'nTracks', 'm2_min', 'm2_max', 'm2_nom', 'm2_final', 'B_0_Theta', 'm2_regreg']#, 'w_min', 'w_max', 'w_nom']
    branch_list +=[ 'b0_TRUEID', 'mu_MC_MOTHER_ID', 'dst_TRUEID', 'b0_TrueHadron_D0_ID', 'd0_M', 'dst_M', 'ctl_reco_regreg', 'ctl_reco_nom', 'ctl_reco_final', 'ctl_reco_max', 'ctl_reco_min', 'ctd_reco_nom', 'ctd_reco_final', 'ctd_reco_max', 'ctd_reco_min','chi_reco_nom', 'chi_reco_final', 'chi_reco_max', 'chi_reco_min', 'w_reco_nom', 'w_reco_final', 'w_reco_max', 'w_reco_min']
    branch_list += mc_weight_base
    branch_list += ['w_true','ctl_true', 'ctd_true', 'chi_true']
    if sample in weights_hammer.keys():
        branch_list += [weights_hammer[sample]]
    df = pd.DataFrame(df_sel.AsNumpy(branch_list))
    df_full = len(df)
    print(df_full)

    df = df.query(tm_strings[sample])
    df['delta_M'] = df['dst_M'] - df['d0_M']
    df = df.query('delta_M >= 140 and delta_M <= 154')


    if sample in weights_hammer.keys():
        df['w_tot'] = df[mc_weight_base + [weights_hammer[sample]]].prod(axis=1)
    else:
        df['w_tot'] = df[mc_weight_base].prod(axis=1)
    print(len(df)/df_full)

    mmiss = 'm2_final'
    yields_pol = []
    errs_pol = []
    if dim == '1D':
        for i in range(len(fit_setup.w_edges[n_w_bins])-1):

            df_w = df.query(f'w_true >= {fit_setup.w_edges[n_w_bins][i]} and w_true < {fit_setup.w_edges[n_w_bins][i+1]}')# and {mmiss} >= {fit_setup.mmiss2_bin_edges[edges][i][0]} and {mmiss} <= {fit_setup.mmiss2_bin_edges[edges][i][-1]}')
            h, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges[edges][i], weights = df_w['w_tot'])
            herr = np.sqrt( np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges[edges][i], weights=df_w['w_tot']**2, density = False)[0])
            yields_pol.append(np.sum(h))
            del df_w
    
    elif dim == '5D':

        bin_no = 0

        for iv0, (bv0l, bv0h) in enumerate(zip(bins[v0], bins[v0][1:])):
            for iv1, (bv1l, bv1h )in enumerate(zip(bins[v1], bins[v1][1:])):
                for bv2l, bv2h in zip(bins[v2], bins[v2][1:]):
                    for bv3l, bv3h in zip(bins[v3], bins[v3][1:]):
                        v0_ = v0.replace('reco_regreg','true')
                        v1_ = v1.replace('reco_regreg','true')
                        v2_ = v2.replace('reco_regreg','true')
                        v3_ = v3.replace('reco_regreg','true')

                        bin_sel = f'{v0_} >= {bv0l} and {v0_} < {bv0h} '
                        bin_sel += f'and {v1_} >= {bv1l} and {v1_} < {bv1h} '
                        bin_sel += f'and {v2_} >= {bv2l} and {v2_} < {bv2h} '
                        bin_sel += f'and {v3_} >= {bv3l} and {v3_} < {bv3h}'
                        mmiss_bins = fit_setup.mmiss2_bin_edges['5D_5'][bin_no]
                        #mmiss_bins = bins['m2miss_reco_approx'][f'{iv0}_{iv1}']
                        #mmiss_bins = fit_setup.mmiss2_bin_edges['1D_5'][iv0]
                        df_w = df.query(bin_sel)
                        h, _ = np.histogram(df_w[mmiss], bins = mmiss_bins, weights = df_w['w_tot'])
                        err = np.sqrt(np.sum(df_w['w_tot']**2))
                        yields_pol.append(np.sum(h))
                        errs_pol.append(err)
                        del df_w
                        bin_no += 1
    yields_up_down[pol] = yields_pol
    errs_up_down[pol] = errs_pol
yields = np.array(yields_up_down['up']) + np.array(yields_up_down['down'])
errs = np.sqrt(np.array(errs_up_down['up'])**2 + np.array(errs_up_down['down'])**2 )/np.sum(yields)
yields /=np.sum(yields)
errs_gen /= np.sum(yields_gen)
yields_gen = np.array(yields_gen)/np.sum(yields_gen)

eff_4d = bfac * yields/yields_gen
err_eff = bfac * eff_4d * np.sqrt((errs_gen/yields_gen)**2 + (errs/yields)**2)


err_dict = {'vals': list(eff_4d),
            'errs' : list(err_eff)}
with open(f"eff_dict_{year}.json", 'w') as f:
    json.dump(err_dict, f, sort_keys=True, indent=4)