"""
Aglomerates some BR values, mostly from PDG, (also has HFLAV https://arxiv.org/pdf/2411.18639)
There are also some theory predictions in https://journals.aps.org/prd/pdf/10.1103/PhysRevD.111.076006

Note that for light feeddown B->D**->D*pi it is predominantly charged B modes that have some D*+ in the chain, as
B+ -> D**0 mu nu and D**0 -> D*+ pi-
By comparison the B0 modes have charged D**+ which often go D**+ -> D*0 pi+

For Ds decays there is:
LHCb paper: https://inspirehep.net/literature/1414195 / https://arxiv.org/pdf/1601.01495
Belle: https://arxiv.org/pdf/2308.08900

Decfile particle names are:
https://gitlab.cern.ch/lhcb-datapkg/Gen/DecFiles/-/blob/v32r39/cmt/decparser/particle_dictionary
"""

import numpy as np
import uncertainties as u


BR = {
    "BR:BdToDstMu"                    : u.ufloat(4.87e-2, 0.09e-2),                     # PDG
    "BR:BdToDstTau"                   : u.ufloat(0.0148, 0.0009),                       # PDG
    "BR:TauToMuNuNu"                  : u.ufloat(0.1739, 0.04),                         # PDG
    "BR:BToD1Lep_D1ToDstPi_HFLAV"     : u.ufloat(0.270e-2, 0.015e-2),                   # HFLAV 2024
    "BR:BToD2Lep_D2ToDstPi_HFLAV"     : u.ufloat(0.085e-2, 0.006e-2),                   # HFLAV 2024
    "BR:BToD1prLep_D1prToDstPi_HFLAV" : u.ufloat(0.185e-2, 0.029e-2),                   # HFLAV 2024
    "BR:BuToD1Lep_D1ToDstpPi"         : u.ufloat(2.84e-3, 0.17e-3),                     # PDG
    "BR:BuToD2Lep_D2ToDstpPi"         : u.ufloat(1.06e-3, 0.18e-3),                     # PDG
    "BR:BuToD1prLep_D1prToDstpPi"     : u.ufloat(1.7e-3, 0.6e-3),                       # PDG
    "BR:BuToDstpPiLep"                : u.ufloat(5.42e-3, 0.28e-3),                     # PDG, inclusive
    "BR:BuToDst0PiPiLep"              : u.ufloat(7.0e-4, 1.7e-4),                       # PDG, inclusive - Babar https://inspirehep.net/literature/1385752
    "BR:BdToD1Lep_D1ToDst0Pi"         : u.ufloat(2.85e-3, 0.25e-3),                     # PDG
    "BR:BdToD2Lep_D2ToDst0Pi"         : u.ufloat(6.6e-4, 1.1e-4),                       # PDG
    "BR:BdToD1prLep_D1prToDst0Pi"     : u.ufloat(2.5e-3, 0.6e-3),                       # PDG
    "BR:BdToDst0PiLep"                : u.ufloat(5.44e-3, 0.28e-3),                     # PDG, inclusive
    "BR:BdToDstpPiPiLep"              : u.ufloat(5.1e-4, 2.3e-4),                       # PDG, inclusive - Belle https://inspirehep.net/literature/2512112
    "BR:BsToDs1Mu_Ds1ToDstKs"         : u.ufloat(2.7e-3, 0.7e-3),                       # PDG, Ds1(2536) - D'_s1 in decfiles. This is Ks only - multiply by 2 for Ks + KL
    "BR:BsToDs2XMu_Ds2ToD0Kp"         : u.ufloat(2.7e-3, 1.0e-3),                       # PDG
    "BR:Ds2ToD0Kp"                    : u.ufloat(0.37, 0.06),                           # PDG
    "R:DststToDstpPi/DststToDst0Pi"   : 0.5,                                            # Isospin symmetry factor
    "R:Ds2ToDstKs/DstToDpKs"          : u.ufloat(0.044, np.sqrt(0.005**2 + 0.011**2)),  # https://arxiv.org/pdf/1601.01495
    "R:Ds2ToDpKs/DstToD0Kp"           : u.ufloat(0.49, np.sqrt(0.1**2 + 0.02**2)),      # https://arxiv.org/pdf/2308.08900
    "R:fs/fd"                         : u.ufloat(0.2539, 0.0079),                        # https://cds.cern.ch/record/2754316/files/LHCb-PAPER-2020-046.pdf?version=2
    "R:BR:BdToDstTau_TauToMuNuNu/BdToDstMu" : u.ufloat(0.281, 0.013), # RD* HFLAV
}

# Some calculations:
BR["BR:BdToDstTau_TauToMuNuNu"] = BR["BR:BdToDstTau"]*BR["BR:TauToMuNuNu"]

BR["BR:Ds2ToDstKs"] = BR["R:Ds2ToDstKs/DstToDpKs"]*BR["R:Ds2ToDpKs/DstToD0Kp"]*BR["BR:Ds2ToD0Kp"]
BR["BR:BsToDs2XMu_Ds2ToDstKs"] = BR["R:Ds2ToDstKs/DstToDpKs"]*BR["R:Ds2ToDpKs/DstToD0Kp"]*BR["BR:BsToDs2XMu_Ds2ToD0Kp"]

BR["BR:BdToD1Lep_D1ToDstpPi"] = BR["R:DststToDstpPi/DststToDst0Pi"]*BR["BR:BdToD1Lep_D1ToDst0Pi"]
BR["BR:BdToD1prLep_D1prToDstpPi"] = BR["R:DststToDstpPi/DststToDst0Pi"]*BR["BR:BdToD1prLep_D1prToDst0Pi"]
BR["BR:BdToD2Lep_D2ToDstpPi"] = BR["R:DststToDstpPi/DststToDst0Pi"]*BR["BR:BdToD2Lep_D2ToDst0Pi"]

BR["BR:BuToDststSum"] = BR["BR:BuToD1Lep_D1ToDstpPi"] + BR["BR:BuToD2Lep_D2ToDstpPi"] + BR["BR:BuToD1prLep_D1prToDstpPi"]
BR["BR:BdToDststSum"] = BR["BR:BdToD2Lep_D2ToDstpPi"] + BR["BR:BdToD1prLep_D1prToDstpPi"] + BR["BR:BdToD1Lep_D1ToDstpPi"]
BR["BR:BsToDsststSum"] = 2*BR["BR:BsToDs1Mu_Ds1ToDstKs"] + 2*BR["BR:BsToDs2XMu_Ds2ToDstKs"]

# We want the ratio of Bs/(B0 + B+). We assume fs/fu = fs/fd ->  fs/(fd+fu) = 1/2*fs/fd
BR["R_BF:BsToDsstst/BToDstst"] = BR["BR:BsToDsststSum"]/(BR["BR:BuToDststSum"] + BR["BR:BdToDststSum"])*(BR["R:fs/fd"])*0.5

# Ratio with considering efficiencies
# Generated number of events sim09k year 2016 magUp + magDown+2017+2018 magUp + magDown
generated_yields = {
    "BdToDstMu"        : 1117276536 + 2226918821,
    "BdToDstTau_TauToMuNuNu": 1263176906 + 1218229494,
    "BuToDststSum"     : 1227691430 + 1175883049+(1140146548+1145355006)+(1301178439+1145355006),
    "BdToDststSum"     : 1199167539 +  1236752126+(1138905547+1131098996)+(1341232837+1131098996),
    "BsToDsststSum"    : 1555761879 + 1443296327 +(1359828427+1368169830)+(1186885466+1368169830)
}

#full q2 range
#template_path = "/eos/lhcb/user/h/havva/q2Measurement-BtoDstMuNu/MCsamples/TrackerOnlyNTuples_Maryland/Angular_RDst_muonic_Run2/tuple_my_variables/weighted_histograms_pid+l0+hlt1_20bins/2016"
template_path = "/eos/lhcb/user/h/havva/q2Measurement-BtoDstMuNu/FitSetup/TotalFit/templates_fullq2_v2"

selected_yields = {
    "BdToDstMu"        : sum(np.loadtxt(f"{template_path}/b0_dstmunu.dat")),
    "BdToDstTau_TauToMuNuNu": sum(np.loadtxt(f"{template_path}/b0_dsttaunu.dat")),
    "BuToDststSum"     : sum(np.loadtxt(f"{template_path}/b0_dststmu_10413.dat"))+sum(np.loadtxt(f"{template_path}/b0_dststmu_20413.dat"))+sum(np.loadtxt(f"{template_path}/b0_dststmu_415.dat")),
    "BdToDststSum"     : sum(np.loadtxt(f"{template_path}/b+_dststmu_10423.dat"))+sum(np.loadtxt(f"{template_path}/b+_dststmu_20423.dat"))+sum(np.loadtxt(f"{template_path}/b+_dststmu_425.dat")),
    "BsToDsststSum"    : sum(np.loadtxt(f"{template_path}/bs_dststmu_10433.dat"))+sum(np.loadtxt(f"{template_path}/bs_dststmu_435.dat"))
}
efficiencies = {}
for key in generated_yields:
    # Here efficiency calculated only for 2016 samples but we use the result as constraint in full Run2 fit
    ## So if we asume similar efficiency for 2017 and 2018 we are fine 
    efficiencies[key] = selected_yields[key]/generated_yields[key] 
    print(f"Efficiency for {key}: {efficiencies[key]:.4e}") 

BR["R_eff:Tau/Mu"] = (efficiencies["BdToDstTau_TauToMuNuNu"]/efficiencies["BdToDstMu"])
BR["R_eff:BsToDsstst/BToDstst"] = (efficiencies["BsToDsststSum"]/(efficiencies["BuToDststSum"] + efficiencies["BdToDststSum"]))
#calculate ratio of yields expected in data
BR["N:BsToDsstst/BToDstst"] = BR["BR:BsToDsststSum"]*efficiencies["BsToDsststSum"]/(BR["BR:BuToDststSum"]*efficiencies["BuToDststSum"] + BR["BR:BdToDststSum"]*efficiencies["BdToDststSum"])*(BR["R:fs/fd"])*0.5
# Print some results
print("R_BF:BsToDsstst/BToDstst", BR["R_BF:BsToDsstst/BToDstst"])
print("R_eff:BsToDsstst/BToDstst", BR["R_eff:BsToDsstst/BToDstst"]) 
print("N(BsToDsstst)/N(BToDstst):", BR["N:BsToDsstst/BToDstst"])

print("N(BdToDstTau_TauToMuNuNu)/N(BdToDstMu)", BR["R:BR:BdToDstTau_TauToMuNuNu/BdToDstMu"]*BR["BR:TauToMuNuNu"]*BR["R_eff:Tau/Mu"])

# print("BR:BdToDstTau_TauToMuNuNu", BR["BR:BdToDstTau_TauToMuNuNu"])
# print("R:Tau/Mu", BR["BR:BdToDstTau_TauToMuNuNu"]/BR["BR:BdToDstMu"])
# print("BR:BuToDststSum", BR["BR:BuToDststSum"])
# print("BR:BdToDststSum", BR["BR:BdToDststSum"])
# print("BR:BsToDsststSum", BR["BR:BsToDsststSum"])
# print("R:BsToDsstst/BToDstst", BR["R:BsToDsstst/BToDstst"])
