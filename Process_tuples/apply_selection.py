import ROOT
from ROOT import TFile, TTree, gROOT, AddressOf, RDataFrame
from math import *
import sys, os
import random

# Open raw tuple
#year = 2016, 2017, 2018
year = sys.argv[1]
#magnet = 'down', 'up'
magnet = sys.argv[2]
#comp_name = 'b0_dstmunu', 'b0_dsttau', 'b0_dststtau', 'b+_dststtau', 'b+_dstdx', 'b0_dstdx', 'b0_dststmu_415', 'b0_dststmu_10413', 'b0_dststmu_20413', 'b+_dststmu_425', 'b+_dststmu_10423', 'b+_dststmu_20423', 'bs_dststmu_10433', 'bs_dststmu_435', 'b+_dststmu_heavy', 'b0_dststmu_heavy'
comp_name = sys.argv[3]
base_path = f"/eos/lhcb/user/h/havva/q2Measurement-BtoDstMuNu/MCsamples/TrackerOnlyNTuples_Maryland/Angular_RDst_muonic_Run2/{year}/"

file_map = {
    "b0_dstmunu": "b0_dstmu/b0_dstmunu",
    "b0_dsttau": "b0_dsttau/b0_dsttaunu",
    "b0_dststtau": "b_dststtau/b0_dststtau",
    "b+_dststtau": "b_dststtau/b+_dststtau",
    "b+_dstdx": "DDX/b+_dstdx",
    "b0_dstdx": "DDX/b0_dstdx",
    "b0_dststmu_415": "b_dststmu/b0_dststmu",
    "b0_dststmu_10413": "b_dststmu/b0_dststmu",
    "b0_dststmu_20413": "b_dststmu/b0_dststmu",
    "b+_dststmu_425": "b_dststmu/b+_dststmu",
    "b+_dststmu_10423": "b_dststmu/b+_dststmu",
    "b+_dststmu_20423": "b_dststmu/b+_dststmu",
    "bs_dststmu_10433": "bs_dststmu/bs_dststmu",
    "bs_dststmu_435": "bs_dststmu/bs_dststmu",
    "b+_dststmu_heavy": "b_dststmu_heavy/b+_dststmu_heavy",
    "b0_dststmu_heavy": "b_dststmu_heavy/b0_dststmu_heavy"
}

# Raw tuple as input
file_dir = file_map[comp_name]+f"_{magnet}"
full_path = base_path + file_dir
fin = TFile.Open(full_path + ".root", "READ")
tin = fin.Get("DecayTree")

# Add necessary firend trees
if year=="2016":
    short_year = 16
elif year=="2017":
    short_year = 17
elif year == "2018":
    short_year = 18

pid_base_path = "/eos/lhcb/user/h/havva/q2Measurement-BtoDstMuNu/PIDCalib2/refcalib_results_ParticlesTogether/"
pid_comp_name = file_map[comp_name].split("/")[-1]
pid_file = f"{pid_comp_name}_{magnet}_PID_eff_Turbo{short_year}_{magnet}.root"
full_pid_path = pid_base_path+pid_file
tin.AddFriend("DecayTreeWpid = PIDCalibTree", full_pid_path)
print(f"[PID] using file: {full_pid_path}")

#HLT1
hlt1_file = TFile.Open(full_path + "_withHLT1.root", "READ")
tin.AddFriend("DecayTreeWhlt1 = TupleB0/DecayTree", hlt1_file)
print(f"[HLT1] using file: {full_path}_withHLT1.root")

#L0muon
l0_file = TFile.Open(full_path + "_withL0emulation.root", "READ")
tin.AddFriend("DecayTreeWl0Muon = DecayTree", l0_file)
print(f"[L0Muon] using file: {full_path}_withL0emulation.root")

df0 = ROOT.RDataFrame(tin)

#Kaon
cutk_PT = 'k_PT > 800' 
cutk_P = 'k_P > 3000'
cutk_IPCHI2 = 'k_IPCHI2_OWNPV > 45'
#cutk_DLLK = 'k_PIDK > 4'  

#Pion                                                                                                   
cutpi_PT = 'pi_PT > 1200'
cutpi_P = 'pi_P > 3000'
cutpi_IPCHI2 = 'pi_IPCHI2_OWNPV > 45'
#cutpi_DLLK = 'pi_PIDK < 2'

#D0
cutd0_PT = 'd0_PT > 2000'
cutd0_IP = 'log(d0_IP_OWNPV) > -3.5'
cutd0_IPCHI2 = 'd0_IPCHI2_OWNPV > 9'
cutd0_FDOWNPV = 'd0_FDCHI2_OWNPV > 250'
cutd0_M = 'd0_M < 1890 and d0_M > 1845'

#Dst
cutdst_M = '(dst_M - d0_M) < 147 and (dst_M - d0_M) > 144' 

#Muon
#cutMuon_ismu = 'mu_isMuon == true'
cutMuon_P = 'mu_P > 6000 && mu_P < 100e3'
#cutMuon_DLL = 'mu_PIDmu > 2.0'
cutMuon_Pt = 'mu_PT > 1200'

#B meson
cutB0_M = '(b0_M) < 5200 and (b0_M) > 2500'

#BDT cut                                                                                                                                                                                    
cutBDT = 'b0_ISOLATION_BDT < 0.15'

#L0 Muon emulation cuts
cutL0mu = 'b0_L0MuonDecision_TOS_emul > 0'

cutB0_OWNPV = '((b0_OWNPV_X)*(b0_OWNPV_X) + (b0_ENDVERTEX_Y - b0_OWNPV_Y)*(b0_ENDVERTEX_Y - b0_OWNPV_Y)) < 7'

# For Tracker Only: Emulated HLT1
#cutHLT1 = '(d0_hlt1_trackmva_tos_emu>0 ||  d0_hlt1_twotrackmva_tos_emu>0)'
cutHLT1 = 'd0_hlt1_trackmva_tos_emu>0 ' #decided on 24/09/2025 since not much loose w/o TwoTrack

# For Tracker Only: Emulated PID 
cutPID = 'PIDCalibEff > 0' # to eliminate -999 (NaN)s

# Add D** selection
# Default: no D** cut
cut_DststID = "true"
if comp_name == "bs_dststmu_10433":
    cut_DststID = "abs(b0_TrueHadron_D0_ID)==10433"
elif comp_name == "bs_dststmu_435":
    cut_DststID = "abs(b0_TrueHadron_D0_ID)==435"
elif comp_name == "b0_dststmu_10413":
    cut_DststID = "abs(b0_TrueHadron_D0_ID)==10413"
elif comp_name == "b0_dststmu_20413":
    cut_DststID = "abs(b0_TrueHadron_D0_ID)==20413"
elif comp_name == "b0_dststmu_415":
    cut_DststID = "abs(b0_TrueHadron_D0_ID)==415"
elif comp_name == "b+_dststmu_10423":
    cut_DststID = "abs(b0_TrueHadron_D0_ID)==10423"
elif comp_name == "b+_dststmu_20423":
    cut_DststID = "abs(b0_TrueHadron_D0_ID)==20423"
elif comp_name == "b+_dststmu_425":
    cut_DststID = "abs(b0_TrueHadron_D0_ID)==425"


df1 = df0.Define("pass_selection_K", "(k_PT > 800)")
df2 = df1.Define("pass_selection_Kpi_d0_Dst", cutk_PT +"&&"+ cutk_P +"&&"+ cutk_IPCHI2 +"&&"+ cutpi_PT +"&&"+ cutpi_P +"&&"+ cutpi_IPCHI2 +"&&"+ cutd0_PT +"&&"+ cutd0_IP +"&&"+ cutd0_IPCHI2 +"&&"+ cutd0_FDOWNPV +"&&"+ cutd0_M +"&&"+ cutdst_M)
df3 = df2.Define("pass_selection_muonP", cutMuon_P)
df4 = df3.Define("pass_selection_muonPT", cutMuon_Pt)
df5 = df4.Define("pass_selection_B0M", cutB0_M)
df6 = df5.Define("pass_selection_IsoBDT", cutBDT)
df7 = df6.Define("pass_selection_L0", cutL0mu)
df8 = df7.Define("pass_selection_HLT1", cutHLT1)
df9 = df8.Define("pass_selection_PID", cutPID)
df10 = df9.Define("pass_selection_DststID", cut_DststID)
#df11 = df9.Define("pass_selection", "pass_selection_Kpi_d0_Dst && pass_selection_muonP && pass_selection_muonPT && pass_selection_IsoBDT && pass_selection_L0 && pass_selection_HLT1 && pass_selection_PID")
df11 = df10.Define("pass_selection", "pass_selection_Kpi_d0_Dst && pass_selection_muonP && pass_selection_muonPT && pass_selection_IsoBDT && pass_selection_L0 && pass_selection_HLT1 && pass_selection_PID && pass_selection_DststID")

print("Writing the new file after cuts ...")
column_name_vector = ROOT.std.vector('string')()
column_name_vector.push_back('pass_selection_Kpi_d0_Dst')
column_name_vector.push_back('pass_selection_muonP')
column_name_vector.push_back('pass_selection_muonPT')
column_name_vector.push_back('pass_selection_B0M')
column_name_vector.push_back('pass_selection_IsoBDT')
column_name_vector.push_back('pass_selection_L0')
column_name_vector.push_back('pass_selection_HLT1')
column_name_vector.push_back('pass_selection_PID')
column_name_vector.push_back('pass_selection')

print("Defined sel variables ...")
#df11.Snapshot("DecayTree", full_path + '_selection_flags.root', column_name_vector)
print(f"[WRITE] output: {year}/{comp_name}_{magnet}_selection_flags.root")
print(f"[INFO] D** cut: {cut_DststID}")

df11.Snapshot("DecayTree", f'{year}/' + comp_name + f'_{magnet}_selection_flags.root', column_name_vector)
print("Done!")