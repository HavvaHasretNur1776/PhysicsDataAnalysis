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
sample = str(sys.argv[1])
year = str(sys.argv[2])
pol = str(sys.argv[3])
dim = str(sys.argv[4])
template_path = fit_setup.template_paths[dim]
postfit_bins = 50
vary_shape = 'nbDstst'
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
postfit_branches = [bpf for bpf in fit_setup.postfit_branch_info.keys()]
if year == 'Run2':
    for i in range(len(fit_setup.w_edges[n_w_bins])-1):
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
        hn = h / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges[dim][i]))
        herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges[dim][i]))
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


    # df = df.Define(
    #     "B_0_Dir_Min",
    #     "ROOT::Math::XYZVector(B_0_Dir).Rotate(B_0_DeltaTheta_Min, B_0_DeltaTheta_Direction)"
    # )

    # df = df.Define(
    #     "B_0_Dir_Max",
    #     "ROOT::Math::XYZVector(B_0_Dir).Rotate(B_0_DeltaTheta_Max, B_0_DeltaTheta_Direction)"
    # )
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

    df_sel = df_sel.Define("ctd_reco_nom","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir, Dst, D0)")
    df_sel = df_sel.Define("ctd_reco_min","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir_Min, Dst, D0)")
    df_sel = df_sel.Define("ctd_reco_max","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir_Max, Dst, D0)")
    df_sel = df_sel.Define("ctd_reco_final","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir_Final, Dst, D0)")

    df_sel = df_sel.Define("chi_reco_nom","ComputeChi(B0_reco_regreg_P,B_0_Dir, Dst, D0, mu, spi, b0_ID)")
    df_sel = df_sel.Define("chi_reco_min","ComputeChi(B0_reco_regreg_P,B_0_Dir_Min, Dst, D0, mu, spi,b0_ID)")
    df_sel = df_sel.Define("chi_reco_max","ComputeChi(B0_reco_regreg_P,B_0_Dir_Max, Dst,mu, spi, D0, b0_ID)")
    df_sel = df_sel.Define("chi_reco_final","ComputeChi(B0_reco_regreg_P,B_0_Dir_Final, Dst, D0,mu, spi, b0_ID)")


    df_sel = df_sel.Define("w_reco_nom","Computew(B0_reco_regreg_P,B_0_Dir, Dst)")
    df_sel = df_sel.Define("w_reco_min","Computew(B0_reco_regreg_P,B_0_Dir_Min, Dst)")
    df_sel = df_sel.Define("w_reco_max","Computew(B0_reco_regreg_P,B_0_Dir_Max, Dst)")
    df_sel = df_sel.Define("w_reco_final","Computew(B0_reco_regreg_P,B_0_Dir_Final, Dst)")
    df_sel = (df_sel
    .Define("b0_ETA", "asinh(b0_PZ/b0_PT)")
    .Define("d0_ETA", "asinh(d0_PZ/d0_PT)")
 
)
    
    branch_list = ['w_reco_regreg', 'm2miss_reco_approx', 'm2_min', 'm2_max', 'm2_nom', 'm2_final', 'B_0_Theta', 'm2_regreg']#'nTracks', 'w_min', 'w_max', 'w_nom']
    branch_list +=[ 'b0_TRUEID', 'mu_MC_MOTHER_ID', 'dst_TRUEID', 'b0_TrueHadron_D0_ID', 'd0_M', 'dst_M', 'ctl_reco_regreg', 'ctl_reco_nom', 'ctl_reco_final', 'ctl_reco_max', 'ctl_reco_min', 'ctd_reco_nom', 'ctd_reco_final', 'ctd_reco_max', 'ctd_reco_min','chi_reco_nom', 'chi_reco_final', 'chi_reco_max', 'chi_reco_min', 'w_reco_nom', 'w_reco_final', 'w_reco_max', 'w_reco_min']
    branch_list += mc_weight_base
    branch_list += postfit_branches

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

    if sample == 'nDstMu':


        plt.hist(df['ctl_reco_nom'], bins = 100,  range = [-1,1],color = 'cornflowerblue', alpha = 0.4, label = 'No smearing')
        comparecounts, comparebins, _ = plt.hist(df['ctl_reco_min'], bins = 100, range = [-1,1], histtype='step', color ='crimson')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = 'v', linestyle = '', color = 'crimson', label = '$\Delta\\theta_-$', ms = 2, lw = 0)
        comparecounts, comparebins, _ = plt.hist(df['ctl_reco_max'], bins = 100, range = [-1,1], histtype='step', color ='crimson')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = '^', linestyle = '', color = 'crimson', label = '$\Delta\\theta_+$', ms = 2, lw = 0)
        comparecounts, comparebins, _ = plt.hist(df['ctl_reco_final'], bins = 100, range = [-1,1], histtype='step', color ='darkgreen')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = 'x', linestyle = '', color = 'darkgreen', label = '$\Delta\\theta_\\alpha$', ms = 2, lw = 0)
        # plt.hist(df['m2_final'], bins = 100, range = [-5,6], histtype='step', label = 'final')
        plt.legend(fontsize = 8)
        plt.savefig(f'ctl_{year}_variation.pdf') 
        plt.clf()

        plt.hist(df['ctd_reco_nom'], bins = 100,  range = [-1,1],color = 'cornflowerblue', alpha = 0.4, label = 'No smearing')
        comparecounts, comparebins, _ = plt.hist(df['ctd_reco_min'], bins = 100, range = [-1,1], histtype='step', color ='crimson')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = 'v', linestyle = '', color = 'crimson', label = '$\Delta\\theta_-$', ms = 2, lw = 0)
        comparecounts, comparebins, _ = plt.hist(df['ctd_reco_max'], bins = 100, range = [-1,1], histtype='step', color ='crimson')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = '^', linestyle = '', color = 'crimson', label = '$\Delta\\theta_+$', ms = 2, lw = 0)
        comparecounts, comparebins, _ = plt.hist(df['ctd_reco_final'], bins = 100, range = [-1,1], histtype='step', color ='darkgreen')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = 'x', linestyle = '', color = 'darkgreen', label = '$\Delta\\theta_\\alpha$', ms = 2, lw = 0)
        # plt.hist(df['m2_final'], bins = 100, range = [-5,6], histtype='step', label = 'final')
        plt.legend(fontsize = 8)
        plt.savefig(f'ctd_{year}_variation.pdf') 
        plt.clf()

        plt.hist(df['chi_reco_nom'], bins = 100,  range = [-np.pi,np.pi],color = 'cornflowerblue', alpha = 0.4, label = 'No smearing')
        comparecounts, comparebins, _ = plt.hist(df['chi_reco_min'], bins = 100, range = [-np.pi,np.pi], histtype='step', color ='crimson')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = 'v', linestyle = '', color = 'crimson', label = '$\Delta\\theta_-$', ms = 2, lw = 0)
        comparecounts, comparebins, _ = plt.hist(df['chi_reco_max'], bins = 100, range = [-np.pi,np.pi], histtype='step', color ='crimson')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = '^', linestyle = '', color = 'crimson', label = '$\Delta\\theta_+$', ms = 2, lw = 0)
        comparecounts, comparebins, _ = plt.hist(df['chi_reco_final'], bins = 100, range = [-np.pi,np.pi], histtype='step', color ='darkgreen')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = 'x', linestyle = '', color = 'darkgreen', label = '$\Delta\\theta_\\alpha$', ms = 2, lw = 0)
        # plt.hist(df['m2_final'], bins = 100, range = [-5,6], histtype='step', label = 'final')
        plt.legend(fontsize = 8)
        plt.savefig(f'chi_{year}_variation.pdf') 
        plt.clf()

        plt.hist(df['w_reco_nom'], bins = 100,  range = [1,1.5],color = 'cornflowerblue', alpha = 0.4, label = 'No smearing')
        comparecounts, comparebins, _ = plt.hist(df['w_reco_min'], bins = 100, range = [1,1.5], histtype='step', color ='crimson')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = 'v', linestyle = '', color = 'crimson', label = '$\Delta\\theta_-$', ms = 2, lw = 0)
        comparecounts, comparebins, _ = plt.hist(df['w_reco_max'], bins = 100, range = [1,1.5], histtype='step', color ='crimson')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = '^', linestyle = '', color = 'crimson', label = '$\Delta\\theta_+$', ms = 2, lw = 0)
        comparecounts, comparebins, _ = plt.hist(df['w_reco_final'], bins = 100, range = [1,1.5], histtype='step', color ='darkgreen')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = 'x', linestyle = '', color = 'darkgreen', label = '$\Delta\\theta_\\alpha$', ms = 2, lw = 0)
        # plt.hist(df['m2_final'], bins = 100, range = [-5,6], histtype='step', label = 'final')
        plt.legend(fontsize = 8)
        plt.savefig(f'w_{year}_variation.pdf') 
        plt.clf()


        plt.hist(df['m2_nom'], bins = 100,  range = [-5,6],color = 'cornflowerblue', alpha = 0.4, label = 'No smearing')
        comparecounts, comparebins, _ = plt.hist(df['m2_min'], bins = 100, range = [-5,6], histtype='step', color ='crimson')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = 'v', linestyle = '', color = 'crimson', label = '$\Delta\\theta_-$', ms = 2, lw = 0)
        comparecounts, comparebins, _ = plt.hist(df['m2_max'], bins = 100, range = [-5,6], histtype='step', color ='crimson')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = '^', linestyle = '', color = 'crimson', label = '$\Delta\\theta_+$', ms = 2, lw = 0)
        comparecounts, comparebins, _ = plt.hist(df['m2_final'], bins = 100, range = [-5,6], histtype='step', color ='darkgreen')
        comparecenters = 0.5 * (comparebins[:-1] + comparebins[1:])
        plt.plot(comparecenters, comparecounts, marker = 'x', linestyle = '', color = 'darkgreen', label = '$\Delta\\theta_\\alpha$', ms = 2, lw = 0)
        plt.legend(fontsize = 8)
        plt.savefig(f'm2_{year}_variation.pdf') 
        plt.clf()
 

    for mmm in ['m2_min', 'm2_nom', 'm2_max', 'm2miss_reco_approx']:
        print(mmm + ' ', df[mmm].min(), ' ', df[mmm].max())
    for l in ['min', 'nom', 'max', 'final']:
        mmiss = 'm2_' + l
        w_smeared = 'w_' + l
        if dim == '1D':
            for i in range(len(fit_setup.w_edges[n_w_bins])-1):

                df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[n_w_bins][i]} and w_reco_regreg < {fit_setup.w_edges[n_w_bins][i+1]} and {mmiss} >= {fit_setup.mmiss2_bin_edges[edges][i][0]} and {mmiss} <= {fit_setup.mmiss2_bin_edges[edges][i][-1]}')
                h, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges[edges][i], weights = df_w['w_tot'])
                herr = np.sqrt( np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges[edges][i], weights=df_w['w_tot']**2, density = False)[0])
                hn, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges[edges][i], weights = df_w['w_tot'], density = True)
                herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges[edges][i]))
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
                del df_w
            df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[n_w_bins][0]} and w_reco_regreg < {fit_setup.w_edges[n_w_bins][-1]}')
            h, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges_full, weights = df_w['w_tot'])
            herr = np.sqrt( np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges[edges][i], weights=df_w['w_tot']**2, density = False)[0])
            hn, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges[edges][i], weights = df_w['w_tot'], density = True)
            herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges[edges][i]))
            f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}{mmiss_lab[l]}_.dat', 'w')
            for j in range(len(hn)):
                a = str(hn[j]) + ' '
                a += str(herrn[j])
                f.write(a+"\n")
            f.close()
            f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}{mmiss_lab[l]}raw.dat', 'w')
            for j in range(len(h)):
                a = str(h[j]) + ' '
                a += str(herr[j])
                f.write(a+"\n")
            f.close()
        
            print(np.sum(h))
        elif dim == '5D' and l == 'final':

            bin_no = 0

            for iv0, (bv0l, bv0h) in enumerate(zip(bins[v0], bins[v0][1:])):
                for iv1, (bv1l, bv1h )in enumerate(zip(bins[v1], bins[v1][1:])):
                    for bv2l, bv2h in zip(bins[v2], bins[v2][1:]):
                        for bv3l, bv3h in zip(bins[v3], bins[v3][1:]):
                            v0_ = v0.replace('regreg','final')
                            v1_ = v1.replace('regreg','final')
                            v2_ = v2.replace('regreg','final')
                            v3_ = v3.replace('regreg','final')
                            bin_sel = f'{v0_} >= {bv0l} and {v0_} < {bv0h} '
                            bin_sel += f'and {v1_} >= {bv1l} and {v1_} < {bv1h} '
                            bin_sel += f'and {v2_} >= {bv2l} and {v2_} < {bv2h} '
                            bin_sel += f'and {v3_} >= {bv3l} and {v3_} < {bv3h}'
                            mmiss_bins = fit_setup.mmiss2_bin_edges['5D_5'][bin_no]
                            #mmiss_bins = bins['m2miss_reco_approx'][f'{iv0}_{iv1}']
                            #mmiss_bins = fit_setup.mmiss2_bin_edges['1D_5'][iv0]
                            df_w = df.query(bin_sel)
                            if sample == 'nDstMu' and (bin_no%15)==0:
                                plt.hist(df_w[mmiss], range = [-1,2.5], bins = 100)
                                plt.savefig(f'check_angular_template.pdf') 
                                plt.clf()
                                plt.hist(df_w[mmiss], bins = 100,range = [mmiss_bins[0], mmiss_bins[-1]])
                                plt.savefig(f'check_angular_template_bin_{bin_no}_{iv0}_{iv1}.pdf') 
                                plt.clf()
                            h, _ = np.histogram(df_w[mmiss], bins = mmiss_bins, weights = df_w['w_tot'])
                            herr = np.sqrt( np.histogram(df_w[mmiss], bins = mmiss_bins, weights=df_w['w_tot']**2, density = False)[0])
                            hn, _ = np.histogram(df_w[mmiss], bins = mmiss_bins, weights = df_w['w_tot'], density = True)
                            herrn = herr / (np.sum(h) * np.diff(mmiss_bins))
                            f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}{mmiss_lab[l]}.dat', 'w')
                            for j in range(len(hn)):

                                a = str(hn[j]) + ' '
                                a += str(herrn[j])
                                f.write(a+"\n")
                            f.close()

                            f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}{mmiss_lab[l]}_raw.dat', 'w')
                            for j in range(len(h)):

                                a = str(h[j]) + ' '
                                a += str(herr[j])
                                f.write(a+"\n")
                            f.close()
                            del df_w
                            bin_no += 1
    for l in postfit_branches:
        mmiss = 'm2_final' 
        w_smeared = 'w_final'
        bin_no = 0

        for iv0, (bv0l, bv0h) in enumerate(zip(bins[v0], bins[v0][1:])):
            for iv1, (bv1l, bv1h )in enumerate(zip(bins[v1], bins[v1][1:])):
                for bv2l, bv2h in zip(bins[v2], bins[v2][1:]):
                    for bv3l, bv3h in zip(bins[v3], bins[v3][1:]):
                        v0_ = v0.replace('regreg','final')
                        v1_ = v1.replace('regreg','final')
                        v2_ = v2.replace('regreg','final')
                        v3_ = v3.replace('regreg','final')
                        bin_sel = f'{v0_} >= {bv0l} and {v0_} < {bv0h} '
                        bin_sel += f'and {v1_} >= {bv1l} and {v1_} < {bv1h} '
                        bin_sel += f'and {v2_} >= {bv2l} and {v2_} < {bv2h} '
                        bin_sel += f'and {v3_} >= {bv3l} and {v3_} < {bv3h}'
                        mmiss_bins = fit_setup.mmiss2_bin_edges['5D_5'][bin_no]
                        #bin_sel += f'and m2_final >= {mmiss_bins[0]} and m2_final < {mmiss_bins[-1]}'
                        bin_sel += f'and m2_final < 0.5'

                        #mmiss_bins = bins['m2miss_reco_approx'][f'{iv0}_{iv1}']
                        #mmiss_bins = fit_setup.mmiss2_bin_edges['1D_5'][iv0]
                        df_w = df.query(bin_sel)

                        h, _ = np.histogram(df_w[l], bins = fit_setup.postfit_branch_info[l][0], range = fit_setup.postfit_branch_info[l][1], weights = df_w['w_tot'])
                        herr = np.sqrt( np.histogram(df_w[mmiss], bins = fit_setup.postfit_branch_info[l][0], range = fit_setup.postfit_branch_info[l][1], weights=df_w['w_tot']**2, density = False)[0])
                        hn, _ = np.histogram(df_w[mmiss], bins = fit_setup.postfit_branch_info[l][0], range = fit_setup.postfit_branch_info[l][1], weights = df_w['w_tot'], density = True)
                        herrn = herr / (np.sum(h) * np.diff(np.linspace(fit_setup.postfit_branch_info[l][1][0], fit_setup.postfit_branch_info[l][1][-1],fit_setup.postfit_branch_info[l][0]+1)))
                        f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}_{l}.dat', 'w')
                        for j in range(len(hn)):

                            a = str(hn[j]) + ' '
                            a += str(herrn[j])
                            f.write(a+"\n")
                        f.close()

                        f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}_{l}_raw.dat', 'w')
                        for j in range(len(h)):

                            a = str(h[j]) + ' '
                            a += str(herr[j])
                            f.write(a+"\n")
                        f.close()
                        del df_w
                        bin_no += 1


elif sample in tuples['data'].keys():
    if sample in ['data','nWSMu']:
        tuple_path = f'/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/sWeighted_data/{year}/{tuples["data"][sample]}_{pol}_selected_withangles_new.root'
        f_main = r.TFile.Open(tuple_path)
        t_main = f_main.Get("DecayTree")
        tuple_friend4 = f'/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/sWeighted_data/{year}/{tuples["data"][sample]}_{pol}_selected_sw_withangles_new.root'
        f_friend4 = r.TFile.Open(tuple_friend4)
        t_friend4 = f_friend4.Get("DecayTree")
        t_main.AddFriend(t_friend4, "sel4")
    else:
        tuple_path = fit_setup.data_tuple_path + f'{year}/{tuples["data"][sample]}_{pol}.root'
        tuple_friend = fit_setup.data_tuple_path + f'{year}/{tuples["data"][sample]}_{pol}_selection_flags.root'
        tuple_friend2 = fit_setup.data_tuple_path + f'{year}/{tuples["data"][sample]}_{pol}_addvars.root'
        tuple_friend3 = fit_setup.data_tuple_path + f'{year}/{tuples["data"][sample]}_{pol}_angles.root'

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
    #df_sel = df_in.Filter("pass_B0_Kinematic == 1 && pass_D0_Kinematic ==1 && pass_Dst_Kinematic ==1 && pass_K_Kinematic ==1 && pass_Mu_Kinematic ==1 && pass_Pi_Kinematic ==1 && pass_HLT1 ==1 && pass_HLT2 ==1 && pass_nSPD_cut == 1 && pass_L0_Had && q2GeV_reco_regreg > 0 and pass_K_PID && pass_Mu_PID && pass_L0_Had && pass_Pi_PID && pass_B0_Isolation == 1")
    if sample in ['data', 'nWSMu']:
        df_sel = df_in
        df_sel = (df_sel
    .Define("b0_ETA", "asinh(b0_PZ/b0_PT)")
    .Define("d0_ETA", "asinh(d0_PZ/d0_PT)")
    .Define("k_ETA", "asinh(k_PZ/k_PT)")
    .Define("pi_ETA", "asinh(pi_PZ/pi_PT)")
    .Define("mu_ETA", "asinh(mu_PZ/mu_PT)")
)
    else:
        df_sel = df_in.Filter("pass_B0_Kinematic == 1 && pass_K_Kinematic ==1 && pass_Mu_Kinematic ==1 && pass_Pi_Kinematic ==1 && pass_HLT1 ==1 && pass_HLT2 ==1 && pass_nSPD_cut == 1 && pass_L0_Had && q2GeV_reco_regreg > 0 and pass_K_PID && pass_Mu_PID && pass_L0_Had && pass_Pi_PID && pass_B0_Isolation == 1")
        df_sel = (df_sel
    .Define("b0_ETA", "asinh(b0_PZ/b0_PT)")
    .Define("d0_ETA", "asinh(d0_PZ/d0_PT)")
    .Define("k_ETA", "asinh(k_PZ/k_PT)")
    .Define("pi_ETA", "asinh(pi_PZ/pi_PT)")
    .Define("mu_ETA", "asinh(mu_PZ/mu_PT)")
)
    
    if sample not in ['data', 'nWSMu']:
        branch_list = ['w_reco_regreg', 'm2miss_reco_approx', 'd0_M', 'dst_M', 'ctl_reco_regreg', 'ctd_reco_regreg', 'chi_reco_regreg']
        branch_list += postfit_branches
        df = pd.DataFrame(df_sel.AsNumpy(branch_list))
        df['delta_M'] = df['dst_M'] - df['d0_M']
        df = df.query('delta_M >= 140 and delta_M <= 154')
    else:
        branch_list = ['w_reco_regreg', 'm2miss_reco_approx','ns_sw', 'ctl_reco_regreg', 'ctd_reco_regreg', 'chi_reco_regreg']
        branch_list += postfit_branches
        df = pd.DataFrame(df_sel.AsNumpy(branch_list))
    
    plt.hist(-df['ctl_reco_regreg'], bins = 100, range = [-1,1], histtype='step', label = 'data')
    plt.legend(fontsize = 8)
    plt.savefig(f'ctl_{year}_min_data.pdf') 
    plt.clf()
    df['ctl_reco_regreg'] = -df['ctl_reco_regreg']
    if dim == '1D':
        for i in range(len(fit_setup.w_edges[n_w_bins])-1):
            df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[n_w_bins][i]} and w_reco_regreg < {fit_setup.w_edges[n_w_bins][i+1]} and m2miss_reco_approx >= {fit_setup.mmiss2_bin_edges[edges][i][0]} and m2miss_reco_approx <= {fit_setup.mmiss2_bin_edges[edges][i][-1]}')
            if sample in ['data', 'nWSMu']:
                h, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges[edges][i], weights = df_w['ns_sw'])
                hnw, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges[edges][i])
                print('Weighted: ', np.sum(h))
                print('no weights: ', np.sum(hnw))
                herr = np.sqrt(h)
                hn, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges[edges][i], density = True, weights = df_w['ns_sw'])
                herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges[edges][i]))
            else:
                h, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges[edges][i])
                herr = np.sqrt(h)
                hn, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges[edges][i], density = True)
                herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges[edges][i]))
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
        df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[n_w_bins][0]} and w_reco_regreg < {fit_setup.w_edges[n_w_bins][-1]}')
        if sample == 'data':
            
            h, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges_full)
            herr = np.sqrt(h)
            hn, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges[edges][i], density = True)
            herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges[edges][i]))
        else:
            h, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges_full)
            herr = np.sqrt(h)
            hn, _ = np.histogram(df_w['m2miss_reco_approx'], bins = fit_setup.mmiss2_bin_edges[edges][i], density = True)
            herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges[edges][i]))
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
    elif dim == '5D':
        if sample == 'data':
            bin_no = 0

            for iv0, (bv0l, bv0h) in enumerate(zip(bins[v0], bins[v0][1:])):
                for iv1, (bv1l, bv1h )in enumerate(zip(bins[v1], bins[v1][1:])):
                    for bv2l, bv2h in zip(bins[v2], bins[v2][1:]):
                        for bv3l, bv3h in zip(bins[v3], bins[v3][1:]):
                            bin_sel = f'{v0} >= {bv0l} and {v0} < {bv0h} '
                            bin_sel += f'and {v1} >= {bv1l} and {v1} < {bv1h} '
                            bin_sel += f'and {v2} >= {bv2l} and {v2} < {bv2h} '
                            bin_sel += f'and {v3} >= {bv3l} and {v3} < {bv3h}'
                            #mmiss_bins = bins['m2miss_reco_approx'][f'{iv0}_{iv1}']
                            mmiss_bins = fit_setup.mmiss2_bin_edges['5D_5'][bin_no]

                            df_w = df.query(bin_sel)
                            h, _ = np.histogram(df_w['m2miss_reco_approx'], bins = mmiss_bins, weights = df_w['ns_sw'])
                            herr = np.sqrt(h)
                            hn, _ = np.histogram(df_w['m2miss_reco_approx'], bins = mmiss_bins, density = True, weights = df_w['ns_sw'])
                            herrn = herr / (np.sum(h) * np.diff(mmiss_bins))
                            f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}.dat', 'w')
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
                            f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}_raw.dat', 'w')
                            for j in range(len(h)):
                                a = str(h[j]) + ' '
                                a += str(herr[j])
                                f.write(a+"\n")
                            f.close()
                            del df_w
                            bin_no +=1
            for l in postfit_branches:
                bin_no = 0

                for iv0, (bv0l, bv0h) in enumerate(zip(bins[v0], bins[v0][1:])):
                    for iv1, (bv1l, bv1h )in enumerate(zip(bins[v1], bins[v1][1:])):
                        for bv2l, bv2h in zip(bins[v2], bins[v2][1:]):
                            for bv3l, bv3h in zip(bins[v3], bins[v3][1:]):
                                bin_sel = f'{v0} >= {bv0l} and {v0} < {bv0h} '
                                bin_sel += f'and {v1} >= {bv1l} and {v1} < {bv1h} '
                                bin_sel += f'and {v2} >= {bv2l} and {v2} < {bv2h} '
                                bin_sel += f'and {v3} >= {bv3l} and {v3} < {bv3h}'
                                #mmiss_bins = bins['m2miss_reco_approx'][f'{iv0}_{iv1}']
                                mmiss_bins = fit_setup.mmiss2_bin_edges['5D_5'][bin_no]
                                #bin_sel += f'and m2miss_reco_approx >= {mmiss_bins[0]} and m2miss_reco_approx < {mmiss_bins[-1]}'
                                bin_sel += f'and m2miss_reco_approx < 0.5'


                                df_w = df.query(bin_sel)
                                h, _ = np.histogram(df_w[l],bins = fit_setup.postfit_branch_info[l][0], range = fit_setup.postfit_branch_info[l][1], weights = df_w['ns_sw'])
                                herr = np.sqrt(h)
                                hn, _ = np.histogram(df_w[l], bins = fit_setup.postfit_branch_info[l][0], range = fit_setup.postfit_branch_info[l][1], density = True, weights = df_w['ns_sw'])
                                herrn = herr / (np.sum(h) * np.diff(np.linspace(fit_setup.postfit_branch_info[l][1][0], fit_setup.postfit_branch_info[l][1][-1],fit_setup.postfit_branch_info[l][0]+1)))
                                f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}_{l}.dat', 'w')
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
                                f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}_{l}_raw.dat', 'w')
                                for j in range(len(h)):
                                    a = str(h[j]) + ' '
                                    a += str(herr[j])
                                    f.write(a+"\n")
                                f.close()
                                del df_w
                                bin_no +=1

        else:
            bin_no = 0
            WSMu_vars = [v0, v1, v2, v3]
            for WSMu_v in WSMu_collapsed_vars:
                WSMu_vars.remove(WSMu_v)
            print('WSMu binning variables: ', WSMu_vars)
            v0 = WSMu_vars[0]
            v1 = WSMu_vars[1]
            v2 = WSMu_collapsed_vars[0]
            v3 = WSMu_collapsed_vars[1]
            for iv0, (bv0l, bv0h) in enumerate(zip(bins[v0], bins[v0][1:])):
                for iv1, (bv1l, bv1h )in enumerate(zip(bins[v1], bins[v1][1:])):
                    for bv2l, bv2h in zip(bins[v2], bins[v2][1:]):
                        for bv3l, bv3h in zip(bins[v3], bins[v3][1:]):
                            bin_sel = f'{v0} >= {bv0l} and {v0} < {bv0h} '
                            bin_sel += f'and {v1} >= {bv1l} and {v1} < {bv1h} '
                            bin_sel += f'and {v2} >= {bins[v2][0]} and {v2} < {bins[v2][-1]} '
                            bin_sel += f'and {v3} >= {bins[v3][0]} and {v3} < {bins[v3][-1]}'
                            #mmiss_bins = bins['m2miss_reco_approx'][f'{iv0}_{iv1}']
                            mmiss_bins = fit_setup.mmiss2_bin_edges['5D_5'][bin_no]
                            
                            #print(bin_sel, ' ', bin_no)
                            df_w = df.query(bin_sel)
                            h, _ = np.histogram(df_w['m2miss_reco_approx'], bins = mmiss_bins, weights =df_w['ns_sw'])
                            herr = np.sqrt(h)
                            hn, _ = np.histogram(df_w['m2miss_reco_approx'], bins = mmiss_bins, density = True, weights = df_w['ns_sw'])
                            herrn = herr / (np.sum(h) * np.diff(mmiss_bins))
                            f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}.dat', 'w')
                            for j in range(len(hn)):
                                a = str(hn[j]) + ' '
                                a += str(herrn[j])
                                f.write(a+"\n")
                            f.close()

                            f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}_raw.dat', 'w')
                            for j in range(len(h)):
                                a = str(h[j]) + ' '
                                a += str(herr[j])
                                f.write(a+"\n")
                            f.close()
                            del df_w
                            bin_no +=1
            for l in postfit_branches:
                bin_no = 0
                WSMu_vars = [v0, v1, v2, v3]
                for WSMu_v in WSMu_collapsed_vars:
                    WSMu_vars.remove(WSMu_v)
                print('WSMu binning variables: ', WSMu_vars)
                v0 = WSMu_vars[0]
                v1 = WSMu_vars[1]
                v2 = WSMu_collapsed_vars[0]
                v3 = WSMu_collapsed_vars[1]
                for iv0, (bv0l, bv0h) in enumerate(zip(bins[v0], bins[v0][1:])):
                    for iv1, (bv1l, bv1h )in enumerate(zip(bins[v1], bins[v1][1:])):
                        for bv2l, bv2h in zip(bins[v2], bins[v2][1:]):
                            for bv3l, bv3h in zip(bins[v3], bins[v3][1:]):
                                bin_sel = f'{v0} >= {bv0l} and {v0} < {bv0h} '
                                bin_sel += f'and {v1} >= {bv1l} and {v1} < {bv1h} '
                                bin_sel += f'and {v2} >= {bins[v2][0]} and {v2} < {bins[v2][-1]} '
                                bin_sel += f'and {v3} >= {bins[v3][0]} and {v3} < {bins[v3][-1]}'
                                #mmiss_bins = bins['m2miss_reco_approx'][f'{iv0}_{iv1}']
                                mmiss_bins = fit_setup.mmiss2_bin_edges['5D_5'][bin_no]
                                #bin_sel += f'and m2miss_reco_approx >= {mmiss_bins[0]} and m2miss_reco_approx < {mmiss_bins[-1]}'
                                bin_sel += f' and m2miss_reco_approx < 0.5'

                                
                                #print(bin_sel, ' ', bin_no)
                                df_w = df.query(bin_sel)
                                h, _ = np.histogram(df_w[l], bins = fit_setup.postfit_branch_info[l][0], range = fit_setup.postfit_branch_info[l][1], weights =df_w['ns_sw'])
                                herr = np.sqrt(h)
                                hn, _ = np.histogram(df_w[l], bins = fit_setup.postfit_branch_info[l][0], range = fit_setup.postfit_branch_info[l][1], density = True, weights = df_w['ns_sw'])
                                herrn = herr / (np.sum(h) * np.diff(np.linspace(fit_setup.postfit_branch_info[l][1][0], fit_setup.postfit_branch_info[l][1][-1],fit_setup.postfit_branch_info[l][0]+1)))
                                f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}_{l}.dat', 'w')
                                for j in range(len(hn)):
                                    a = str(hn[j]) + ' '
                                    a += str(herrn[j])
                                    f.write(a+"\n")
                                f.close()

                                f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{bin_no}_{l}_raw.dat', 'w')
                                for j in range(len(h)):
                                    a = str(h[j]) + ' '
                                    a += str(herr[j])
                                    f.write(a+"\n")
                                f.close()
                                del df_w
                                bin_no +=1




else:
    subtemplates = {}
    for subsample in tm_strings[sample]:
        subtemplates[subsample] = {}
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
        tuple_friend10 = fit_setup.tuple_path + f'{year}/{tuples["MC"][sample][subsample]}_{pol}_pidcalib.root'
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
        f_friend10 = r.TFile.Open(tuple_friend10)

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
        t_friend10 = f_friend10.Get("DecayTree")
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
        t_main.AddFriend(t_friend10, "sel10")

        df_in = r.RDataFrame(t_main)

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


        # df = df.Define(
        #     "B_0_Dir_Min",
        #     "ROOT::Math::XYZVector(B_0_Dir).Rotate(B_0_DeltaTheta_Min, B_0_DeltaTheta_Direction)"
        # )

        # df = df.Define(
        #     "B_0_Dir_Max",
        #     "ROOT::Math::XYZVector(B_0_Dir).Rotate(B_0_DeltaTheta_Max, B_0_DeltaTheta_Direction)"
        # )
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
        .Define("m2_min", "ComputeM2(B_0_Dir_Min, Dst, D0, mu, B0_reco_approx_P)")
        .Define("m2_max", "ComputeM2(B_0_Dir_Max, Dst, D0, mu, B0_reco_approx_P)")
        .Define("m2_final", "ComputeM2(B_0_Dir_Final, Dst, D0, mu, B0_reco_approx_P)"))
        df_sel = df_sel.Define("ctl_reco_nom","ComputeCosThetaL(B0_reco_regreg_P,B_0_Dir, Dst, mu)")
        df_sel = df_sel.Define("ctl_reco_min","ComputeCosThetaL(B0_reco_regreg_P,B_0_Dir_Min, Dst, mu)")
        df_sel = df_sel.Define("ctl_reco_max","ComputeCosThetaL(B0_reco_regreg_P,B_0_Dir_Max, Dst, mu)")
        df_sel = df_sel.Define("ctl_reco_final","ComputeCosThetaL(B0_reco_regreg_P,B_0_Dir_Final, Dst, mu)")

        df_sel = df_sel.Define("ctd_reco_nom","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir, Dst, D0)")
        df_sel = df_sel.Define("ctd_reco_min","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir_Min, Dst, D0)")
        df_sel = df_sel.Define("ctd_reco_max","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir_Max, Dst, D0)")
        df_sel = df_sel.Define("ctd_reco_final","ComputeCosThetaD(B0_reco_regreg_P,B_0_Dir_Final, Dst, D0)")

        df_sel = df_sel.Define("chi_reco_nom","ComputeChi(B0_reco_regreg_P,B_0_Dir, Dst, D0, mu, spi, b0_ID)")
        df_sel = df_sel.Define("chi_reco_min","ComputeChi(B0_reco_regreg_P,B_0_Dir_Min, Dst, D0, mu, spi,b0_ID)")
        df_sel = df_sel.Define("chi_reco_max","ComputeChi(B0_reco_regreg_P,B_0_Dir_Max, Dst,mu, spi, D0, b0_ID)")
        df_sel = df_sel.Define("chi_reco_final","ComputeChi(B0_reco_regreg_P,B_0_Dir_Final, Dst, D0,mu, spi, b0_ID)")


        df_sel = df_sel.Define("w_reco_nom","Computew(B0_reco_regreg_P,B_0_Dir, Dst)")
        df_sel = df_sel.Define("w_reco_min","Computew(B0_reco_regreg_P,B_0_Dir_Min, Dst)")
        df_sel = df_sel.Define("w_reco_max","Computew(B0_reco_regreg_P,B_0_Dir_Max, Dst)")
        df_sel = df_sel.Define("w_reco_final","Computew(B0_reco_regreg_P,B_0_Dir_Final, Dst)")
        df_sel = (df_sel
    .Define("b0_ETA", "asinh(b0_PZ/b0_PT)")
    .Define("d0_ETA", "asinh(d0_PZ/d0_PT)")

)
       

        branch_list = ['w_reco_regreg', 'm2miss_reco_approx', 'm2_nom', 'm2_min', 'm2_max', 'm2_final']
        branch_list +=[ 'b0_TRUEID', 'mu_MC_MOTHER_ID', 'dst_TRUEID', 'b0_TrueHadron_D0_ID', 'd0_M', 'dst_M','ctl_reco_regreg', 'ctl_reco_nom', 'ctl_reco_final', 'ctl_reco_max', 'ctl_reco_min', 'ctd_reco_nom', 'ctd_reco_final', 'ctd_reco_max', 'ctd_reco_min','chi_reco_nom', 'chi_reco_final', 'chi_reco_max', 'chi_reco_min', 'w_reco_nom', 'w_reco_final', 'w_reco_max', 'w_reco_min']
        branch_list += mc_weight_base
        branch_list += postfit_branches

        if sample in weights_hammer.keys():
            branch_list += [weights_hammer[sample][subsample]]
        df = pd.DataFrame(df_sel.AsNumpy(branch_list))

        df_full = len(df)
        print(df_full)

        df = df.query(tm_strings[sample][subsample])
        df['delta_M'] = df['dst_M'] - df['d0_M']
        df = df.query('delta_M >= 140 and delta_M <= 154')
        #plt.hist(df['ctl_reco_nom'], range = [-1,1], bins = 100, alpha = 0.5)
        plt.hist(df['ctl_reco_final'], range = [-1,1], bins = 100, alpha = 0.5)
        plt.xlabel('$\\cos\\theta_L$')
        plt.savefig(f'regreg_vs_final_{year}_{subsample}_{pol}.pdf') 
        plt.clf()

        xp = df["ctl_reco_final"]
        yp = df["w_reco_final"]

        xbins = [-1, -0.3543,0,1] 
        ybins = [1.  , 1.1 , 1.2,  1.3 , 1.4 ,  1.51 ]
        Hp, xedgesp, yedgesp = np.histogram2d(xp, yp, bins=[xbins, ybins])

    
        fig, ax = plt.subplots(figsize=(8, 6))

        mesh = ax.pcolormesh(xedgesp, yedgesp, Hp.T, cmap="viridis")
        plt.colorbar(mesh, ax=ax, label="Counts")
        
        for i in range(len(xedgesp) - 1):
            for j in range(len(yedgesp) - 1):
                xc = (xedgesp[i] + xedgesp[i+1]) / 2
                yc = (yedgesp[j] + yedgesp[j+1]) / 2
                ax.text(
                    xc, yc,
                    f"{int(Hp[i, j])}",
                    ha="center", va="center",
                    color="black" if Hp[i, j] > Hp.max()/2 else "white"
                )

        ax.set_xlabel("$\\cos\\theta_L$")
        ax.set_ylabel("$w$")

        plt.savefig(f'2Dpop_{year}_{subsample}_{pol}.pdf') 
        plt.clf()

        # print(df.head(10))
        if sample in weights_hammer.keys():
            df['w_tot'] = df[mc_weight_base + [weights_hammer[sample][subsample]]].prod(axis=1)
        else:
            df['w_tot'] = df[mc_weight_base].prod(axis=1)
        print(len(df)/df_full)
        sel_evts = len(df)
        for l in ['min', 'nom', 'max', 'final']:
            mmiss = 'm2_' + l
            w_smeared = 'w_reco' + l
            subtemplates[subsample][l] = {}
            if dim == '1D':
                for i in range(len(fit_setup.w_edges[n_w_bins])-1):
                    df_w = df.query(f'w_reco_regreg >= {fit_setup.w_edges[n_w_bins][i]} and w_reco_regreg < {fit_setup.w_edges[n_w_bins][i+1]}  and {mmiss} >= {fit_setup.mmiss2_bin_edges[edges][i][0]} and {mmiss} <= {fit_setup.mmiss2_bin_edges[edges][i][-1]}')
                    h, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges[edges][i], weights = df_w['w_tot'])
                    hn, _ = np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges[edges][i], weights = df_w['w_tot'], density = True)
                    herr = np.sqrt( np.histogram(df_w[mmiss], bins = fit_setup.mmiss2_bin_edges[edges][i], weights=df_w['w_tot']**2, density = False)[0])
                    herrn = herr / (np.sum(h) * np.diff(fit_setup.mmiss2_bin_edges[edges][i]))
                    subtemplates[subsample][l][i] = {'nev': np.sum(h),'template': hn, 'template_err': herrn, 'contrib': rel_BFs[sample][subsample][0] *(sel_evts) / (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample]), 'contrib_error': np.sqrt((rel_BFs[sample][subsample][0]**2  * sel_evts**2) * ((rel_BFs[sample][subsample][1]/ rel_BFs[sample][subsample][0])**2 + (np.sqrt(sel_evts)/sel_evts)**2) )/ (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample])}
                    del df_w
            elif dim == '5D' and l == 'final':
                bin_no = 0
                df_plot = df.query('ctl_reco_final >= -1 and ctl_reco_final < -0.3543 and w_reco_final >= 1.4 and w_reco_final < 1.51')
                plt.hist(df_plot['w_tot'], bins = 20)
                plt.savefig(f'weights_w4ctl0_{year}_{subsample}_{pol}.pdf') 
                plt.clf()
                for iv0, (bv0l, bv0h) in enumerate(zip(bins[v0], bins[v0][1:])):
                    for iv1, (bv1l, bv1h )in enumerate(zip(bins[v1], bins[v1][1:])):
                        for bv2l, bv2h in zip(bins[v2], bins[v2][1:]):
                            for bv3l, bv3h in zip(bins[v3], bins[v3][1:]):
                                v0_ = v0.replace('regreg','final')
                                v1_ = v1.replace('regreg','final')
                                v2_ = v2.replace('regreg','final')
                                v3_ = v3.replace('regreg','final')
                                bin_sel = f'{v0_} >= {bv0l} and {v0_} < {bv0h} '
                                bin_sel += f'and {v1_} >= {bv1l} and {v1_} < {bv1h} '
                                bin_sel += f'and {v2_} >= {bv2l} and {v2_} < {bv2h} '
                                bin_sel += f'and {v3_} >= {bv3l} and {v3_} < {bv3h}'
                                #mmiss_bins = bins['m2miss_reco_approx'][f'{iv0}_{iv1}']
                                mmiss_bins = fit_setup.mmiss2_bin_edges['5D_5'][bin_no]
                                df_w = df.query(bin_sel)
                                h, _ = np.histogram(df_w[mmiss], bins = mmiss_bins, weights = df_w['w_tot'])
                                hn, _ = np.histogram(df_w[mmiss], bins = mmiss_bins, weights = df_w['w_tot'], density = True)
                                herr = np.sqrt( np.histogram(df_w[mmiss], bins = mmiss_bins, weights=df_w['w_tot']**2, density = False)[0])
                                herrn = herr / (np.sum(h) * np.diff(mmiss_bins))
                                subtemplates[subsample][l][bin_no] = {'nev': np.sum(np.nan_to_num(h)),'template': np.nan_to_num(hn), 'template_err': np.nan_to_num(herrn), 'contrib': rel_BFs[sample][subsample][0] *(sel_evts) / (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample]), 'contrib_error': np.sqrt((rel_BFs[sample][subsample][0]**2  * sel_evts**2) * ((rel_BFs[sample][subsample][1]/ rel_BFs[sample][subsample][0])**2 + (np.sqrt(sel_evts)/sel_evts)**2) )/ (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample])}
                                #subtemplates_raw[subsample][bin_no] = {'nev': np.sum(h),'template': h, 'template_err': herr, 'contrib': rel_BFs[sample][subsample][0] *(sel_evts) / (fit_setup.gen_yields[year][sample][subsample]* dec_fracs[sample][subsample]), 'contrib_error': np.sqrt((rel_BFs[sample][subsample][0]**2  * sel_evts**2) * ((rel_BFs[sample][subsample][1]/ rel_BFs[sample][subsample][0])**2 + (np.sqrt(sel_evts)/sel_evts)**2) )/ (fit_setup.gen_yields[year][sample][subsample]* dec_fracs[sample][subsample])}

                                del df_w
                                bin_no +=1
    #print(subtemplates)
        for l in postfit_branches:
            mmiss = 'm2_final' 
            w_smeared = 'w_reco_final' 
            subtemplates[subsample][l] = {}
            bin_no = 0
            for iv0, (bv0l, bv0h) in enumerate(zip(bins[v0], bins[v0][1:])):
                for iv1, (bv1l, bv1h )in enumerate(zip(bins[v1], bins[v1][1:])):
                    for bv2l, bv2h in zip(bins[v2], bins[v2][1:]):
                        for bv3l, bv3h in zip(bins[v3], bins[v3][1:]):
                            v0_ = v0.replace('regreg','final')
                            v1_ = v1.replace('regreg','final')
                            v2_ = v2.replace('regreg','final')
                            v3_ = v3.replace('regreg','final')
                            bin_sel = f'{v0_} >= {bv0l} and {v0_} < {bv0h} '
                            bin_sel += f'and {v1_} >= {bv1l} and {v1_} < {bv1h} '
                            bin_sel += f'and {v2_} >= {bv2l} and {v2_} < {bv2h} '
                            bin_sel += f'and {v3_} >= {bv3l} and {v3_} < {bv3h}'
                            #mmiss_bins = bins['m2miss_reco_approx'][f'{iv0}_{iv1}']
                            mmiss_bins = fit_setup.mmiss2_bin_edges['5D_5'][bin_no]
                            #bin_sel += f'and {mmiss} >= {mmiss_bins[0]} and {mmiss} < {mmiss_bins[-1]}'
                            bin_sel += f'and {mmiss} < 0.5'
                            df_w = df.query(bin_sel)
                            h, _ = np.histogram(df_w[l],  bins = fit_setup.postfit_branch_info[l][0], range = fit_setup.postfit_branch_info[l][1],  weights = df_w['w_tot'])
                            hn, _ = np.histogram(df_w[l],  bins = fit_setup.postfit_branch_info[l][0], range = fit_setup.postfit_branch_info[l][1],  weights = df_w['w_tot'], density = True)
                            herr = np.sqrt( np.histogram(df_w[l],  bins = fit_setup.postfit_branch_info[l][0], range = fit_setup.postfit_branch_info[l][1],  weights=df_w['w_tot']**2, density = False)[0])
                            herrn = herr / (np.sum(h) * np.diff(np.linspace(fit_setup.postfit_branch_info[l][1][0], fit_setup.postfit_branch_info[l][1][-1],fit_setup.postfit_branch_info[l][0]+1)))
                            subtemplates[subsample][l][bin_no] = {'nev': np.sum(np.nan_to_num(h)),'template': np.nan_to_num(hn), 'template_err': np.nan_to_num(herrn), 'contrib': rel_BFs[sample][subsample][0] *(sel_evts) / (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample]), 'contrib_error': np.sqrt((rel_BFs[sample][subsample][0]**2  * sel_evts**2) * ((rel_BFs[sample][subsample][1]/ rel_BFs[sample][subsample][0])**2 + (np.sqrt(sel_evts)/sel_evts)**2) )/ (fit_setup.gen_yields[year][sample][subsample][pol]* dec_fracs[sample][subsample])}
                            #subtemplates_raw[subsample][bin_no] = {'nev': np.sum(h),'template': h, 'template_err': herr, 'contrib': rel_BFs[sample][subsample][0] *(sel_evts) / (fit_setup.gen_yields[year][sample][subsample]* dec_fracs[sample][subsample]), 'contrib_error': np.sqrt((rel_BFs[sample][subsample][0]**2  * sel_evts**2) * ((rel_BFs[sample][subsample][1]/ rel_BFs[sample][subsample][0])**2 + (np.sqrt(sel_evts)/sel_evts)**2) )/ (fit_setup.gen_yields[year][sample][subsample]* dec_fracs[sample][subsample])}

                            del df_w
                            bin_no +=1
    for l in ['min', 'nom', 'max', 'final']:
        if dim == '1D':
            for i in range(len(fit_setup.w_edges[n_w_bins])-1):
                
                full_temp = np.zeros(len(fit_setup.mmiss2_bin_edges[edges][i])-1)
                full_err = np.zeros(len(fit_setup.mmiss2_bin_edges[edges][i])-1)
                num_err = np.zeros(len(fit_setup.mmiss2_bin_edges[edges][i])-1)
                den_err = np.zeros(len(fit_setup.mmiss2_bin_edges[edges][i])-1)
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
                full_temp /= full_contrib

                f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}_bin{i}{mmiss_lab[l]}.dat', 'w')
                for j in range(len(full_temp)):
                    a = str(full_temp[j]) + ' '
                    a += str(full_err[j])
                    f.write(a+"\n")
                f.close()
                full_temp *= nev * np.diff(fit_setup.mmiss2_bin_edges[edges][i])
                full_err *= nev * np.diff(fit_setup.mmiss2_bin_edges[edges][i])
                f = open(fit_setup.template_path + f'/{year}/' + f'{template_label[sample]}_{pol}_bin{i}{mmiss_lab[l]}_raw.dat', 'w')
                for j in range(len(full_temp)):
                    a = str(full_temp[j]) + ' '
                    a += str(full_err[j])
                    f.write(a+"\n")
                f.close()
        elif dim == '5D' and l == 'final':
            for k in subtemplates.keys():
                print(k, ': ', subtemplates[k][l][0])
            ntemplates = nbins[v0] * nbins[v1] * nbins[v2] * nbins[v3]
            for i in range(ntemplates):
                full_temp = np.zeros(len(fit_setup.mmiss2_bin_edges[edges][i])-1)
                full_err = np.zeros(len(fit_setup.mmiss2_bin_edges[edges][i])-1)
                num_err = np.zeros(len(fit_setup.mmiss2_bin_edges[edges][i])-1)
                den_err = np.zeros(len(fit_setup.mmiss2_bin_edges[edges][i])-1)
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

                f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{i}{mmiss_lab[l]}.dat', 'w')
                for j in range(len(full_temp)):
                    a = str(full_temp[j]) + ' '
                    a += str(full_err[j])
                    f.write(a+"\n")
                f.close()
                full_temp *= nev * np.diff(fit_setup.mmiss2_bin_edges[edges][i])
                full_err *= nev * np.diff(fit_setup.mmiss2_bin_edges[edges][i])
                f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{i}{mmiss_lab[l]}_raw.dat', 'w')
                for j in range(len(full_temp)):
                    a = str(full_temp[j]) + ' '
                    a += str(full_err[j])
                    f.write(a+"\n")
                f.close()


    for l in postfit_branches:
        for k in subtemplates.keys():
            print(k, ': ', subtemplates[k][l][0])
        ntemplates = nbins[v0] * nbins[v1] * nbins[v2] * nbins[v3]
        for i in range(ntemplates):
            full_temp = np.zeros(fit_setup.postfit_branch_info[l][0])
            full_err = np.zeros(fit_setup.postfit_branch_info[l][0])
            num_err = np.zeros(fit_setup.postfit_branch_info[l][0])
            den_err = np.zeros(fit_setup.postfit_branch_info[l][0])
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

            f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{i}_{l}.dat', 'w')
            for j in range(len(full_temp)):
                a = str(full_temp[j]) + ' '
                a += str(full_err[j])
                f.write(a+"\n")
            f.close()
            full_temp *= nev * np.diff(np.linspace(fit_setup.postfit_branch_info[l][1][0], fit_setup.postfit_branch_info[l][1][-1],fit_setup.postfit_branch_info[l][0]+1))
            full_err *= nev * np.diff(np.linspace(fit_setup.postfit_branch_info[l][1][0], fit_setup.postfit_branch_info[l][1][-1],fit_setup.postfit_branch_info[l][0]+1))
            f = open(fit_setup.template_path_5D + f'{year}/' + f'{template_label[sample]}_{pol}_bin{i}_{l}_raw.dat', 'w')
            for j in range(len(full_temp)):
                a = str(full_temp[j]) + ' '
                a += str(full_err[j])
                f.write(a+"\n")
            f.close()