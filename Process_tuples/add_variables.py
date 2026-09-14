import ROOT
from ROOT import TFile, TTree, gROOT, AddressOf, RDataFrame
from math import *
import sys, os
import random

# Open raw tuple
year = sys.argv[1]
magnet = sys.argv[2]
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

# Add selections tree
#_name_selection_flags =  file_map[comp_name]+f"_{magnet}_selection_flags.root"
_name_selection_flags =  f"{year}/" + comp_name + f"_{magnet}_selection_flags.root"
tin.AddFriend("DecayTree_wpasscut = DecayTree", _name_selection_flags)
print(f"[SELECTION FLAGS] friend added using file: {_name_selection_flags}")

# Add necessary friend trees to record pid weights
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
print(f"[PID] friend added using file: {full_pid_path}")

#outfilename = file_map[comp_name] +f"_{magnet}_allVariables.root"
outfilename = f"{year}/" + comp_name +f"_{magnet}_allVariables.root"


auxiliar_functions_code =''' 

double get_costhetal(TLorentzVector b0_PXYZE, TLorentzVector mu_PXYZE, TLorentzVector W_PXYZE)
{
  double costheta_l;
  TVector3 b0_boost = -1.0*(b0_PXYZE.BoostVector());
  mu_PXYZE.Boost(b0_boost);
  W_PXYZE.Boost(b0_boost);
  b0_PXYZE.Boost(b0_boost);
  TVector3 W_boost = -1.0*(W_PXYZE.BoostVector());
  mu_PXYZE.Boost(W_boost);
  b0_PXYZE.Boost(W_boost);

  costheta_l = TMath::Cos((-1.0*b0_PXYZE.Vect().Unit()).Angle(mu_PXYZE.Vect().Unit()));
  return (costheta_l);
}

double get_thetal(TLorentzVector b0_PXYZE, TLorentzVector mu_PXYZE, TLorentzVector W_PXYZE)
{
  double theta_l;
  TVector3 b0_boost = -1.0*(b0_PXYZE.BoostVector());
  mu_PXYZE.Boost(b0_boost);
  W_PXYZE.Boost(b0_boost);
  b0_PXYZE.Boost(b0_boost);
  TVector3 W_boost = -1.0*(W_PXYZE.BoostVector());
  mu_PXYZE.Boost(W_boost);
  b0_PXYZE.Boost(W_boost);

  theta_l = (-1.0*b0_PXYZE.Vect().Unit()).Angle(mu_PXYZE.Vect().Unit());
  return (theta_l);
}
double get_sinthetal(TLorentzVector b0_PXYZE, TLorentzVector mu_PXYZE, TLorentzVector W_PXYZE)
{
  double sintheta_l;
  TVector3 b0_boost = -1.0*(b0_PXYZE.BoostVector());
  mu_PXYZE.Boost(b0_boost);
  W_PXYZE.Boost(b0_boost);
  b0_PXYZE.Boost(b0_boost);
  TVector3 W_boost = -1.0*(W_PXYZE.BoostVector());
  mu_PXYZE.Boost(W_boost);
  b0_PXYZE.Boost(W_boost);

  sintheta_l = TMath::Sin((-1.0*b0_PXYZE.Vect().Unit()).Angle(mu_PXYZE.Vect().Unit()));
  return (sintheta_l);
}


double get_costhetad(TLorentzVector b0_PXYZE, TLorentzVector Dp_PXYZE, TLorentzVector dst_PXYZE, TLorentzVector W_PXYZE)
{
  double costheta_d;
  //boost in D* rest frame                                                                                                                                                                   
  dst_PXYZE.Boost(-b0_PXYZE.BoostVector());
  Dp_PXYZE.Boost(-b0_PXYZE.BoostVector());
  b0_PXYZE.Boost(-b0_PXYZE.BoostVector());
  Dp_PXYZE.Boost(-dst_PXYZE.BoostVector());
  b0_PXYZE.Boost(-dst_PXYZE.BoostVector());

  costheta_d = TMath::Cos(Dp_PXYZE.Vect().Unit().Angle(b0_PXYZE.Vect().Unit()));

  return (costheta_d);
}

double get_chi(TLorentzVector nu_PXYZE, TLorentzVector mu_PXYZE, TLorentzVector pi_PXYZE, TLorentzVector Dp_PXYZE, TLorentzVector b0_PXYZE, TLorentzVector W_PXYZE, TLorentzVector dst_PXYZE)
{
  double phi_p,s,c,phi_p_test,sin,cos,phi_p_test_2;
  W_PXYZE.Boost(-b0_PXYZE.BoostVector());
  dst_PXYZE.Boost(-b0_PXYZE.BoostVector());
  nu_PXYZE.Boost(-b0_PXYZE.BoostVector());
  mu_PXYZE.Boost(-b0_PXYZE.BoostVector());
  pi_PXYZE.Boost(-b0_PXYZE.BoostVector());
  Dp_PXYZE.Boost(-b0_PXYZE.BoostVector());
  b0_PXYZE.Boost(-b0_PXYZE.BoostVector());
  TVector3 vec1(Dp_PXYZE.Vect().Cross(pi_PXYZE.Vect())); //unit vector normal to the plane of D and pi                                                                                      
  TVector3 vec2(mu_PXYZE.Vect().Cross(nu_PXYZE.Vect())); //unit vector normal to the plae of mu and nu                                                                                      
  phi_p = vec2.Unit().Angle(vec1.Unit());
  c = vec2.Unit().Dot(vec1.Unit());
  s = (vec2.Unit().Cross(vec1.Unit())).Mag();
  phi_p_test = TMath::ATan2(s,c);

  TVector3 vec_DstD0(dst_PXYZE.Vect().Unit().Cross(Dp_PXYZE.Vect().Unit()));
  TVector3 vec_WMu(W_PXYZE.Vect().Unit().Cross(mu_PXYZE.Vect().Unit()));
  TVector3 vec_y(vec_DstD0.Unit().Cross(dst_PXYZE.Vect()));

cos = vec_WMu.Dot(vec_DstD0.Unit());
  sin = vec_WMu.Dot(vec_y.Unit());
  phi_p_test_2 = TMath::ATan2(sin,cos);

  return (phi_p_test_2);
}

TLorentzVector get_Bmom_rec_Brestframe(TVector3 DIRECTION, TLorentzVector dst_PXYZE, TLorentzVector mu_PXYZE, double Bplus_M, double Bplus_PZ)
{
  double norm =  TMath::Sqrt(DIRECTION.Mag2());
  double b0_cosangle_X = DIRECTION.X()/norm;
  double b0_cosangle_Y = DIRECTION.Y()/norm;
  double b0_cosangle_Z = DIRECTION.Z()/norm;
  double b0_angle_Z = TMath::ACos(b0_cosangle_Z);
                                                                                                                                                                      
  double b0_M_val = 5279.610;
  double Dstmu_M = TMath::Sqrt((dst_PXYZE.E() + mu_PXYZE.E())*(dst_PXYZE.E() + mu_PXYZE.E()) - (dst_PXYZE.X() + mu_PXYZE.X())*(dst_PXYZE.X() + mu_PXYZE.X()) - (dst_PXYZE.Y() + mu_PXYZE.Y())*(dst_PXYZE.Y() + mu_PXYZE.Y()) - (dst_PXYZE.Z() + mu_PXYZE.Z())*(dst_PXYZE.Z() + mu_PXYZE.Z()));                                                     
  double Dstmu_PZ = dst_PXYZE.Z() + mu_PXYZE.Z();
  double b0_PZ_rec = b0_M_val/Bplus_M * Bplus_PZ;
  double b0_P_rec = b0_PZ_rec * (TMath::Sqrt(1.0 + TMath::Tan(b0_angle_Z) * TMath::Tan(b0_angle_Z)));
  double b0_PX_rec = b0_cosangle_X * b0_P_rec;
  double b0_PY_rec = b0_cosangle_Y * b0_P_rec;
  double b0_PZ_rec_test = b0_cosangle_Z * b0_P_rec;
  double b0_PE_rec = TMath::Sqrt((b0_M_val * b0_M_val) + (b0_P_rec * b0_P_rec));
  TLorentzVector PXYZE_rec(b0_PX_rec, b0_PY_rec, b0_PZ_rec_test, b0_PE_rec);

  return (PXYZE_rec);
}

double get_b0_flight_angle(TVector3 DIRECTION)
{
  double norm =  TMath::Sqrt(DIRECTION.Mag2());
  double b0_cosangle_Z = DIRECTION.Z()/norm;
  double b0_angle_Z = TMath::ACos(b0_cosangle_Z);
  return b0_angle_Z;
}

TLorentzVector get_Bmom_rec_quadratic(TVector3 DIRECTION, TLorentzVector dst_PXYZE, TLorentzVector mu_PXYZE, TLorentzVector Pvis_PXYZE, int sign)
{
  double norm =  TMath::Sqrt(DIRECTION.Mag2());
  double b0_cosangle_X = DIRECTION.X()/norm;
  double b0_cosangle_Y = DIRECTION.Y()/norm;
  double b0_cosangle_Z = DIRECTION.Z()/norm;
  double theta_vis = DIRECTION.Unit().Angle(Pvis_PXYZE.Vect().Unit());    
  double P_long = TMath::Cos(theta_vis) * TMath::Sqrt(Pvis_PXYZE.Vect().Mag2());
  double P_trans = TMath::Sin(theta_vis) * TMath::Sqrt(Pvis_PXYZE.Vect().Mag2());
  double b0_M_val = 5279.610;
  double b0_M_val2 = b0_M_val * b0_M_val;
  double Dstmu_M2 = (dst_PXYZE.E() + mu_PXYZE.E())*(dst_PXYZE.E() + mu_PXYZE.E()) - (dst_PXYZE.X() + mu_PXYZE.X())*(dst_PXYZE.X() + mu_PXYZE.X()) - (dst_PXYZE.Y() + mu_PXYZE.Y())*(dst_PXYZE.Y() + mu_PXYZE.Y()) - (dst_PXYZE.Z() + mu_PXYZE.Z())*(dst_PXYZE.Z() + mu_PXYZE.Z());
  double a1 = b0_M_val2 - Dstmu_M2 - (2*(P_trans * P_trans));
  double a2 = (P_long * P_long) - (Pvis_PXYZE.E() * Pvis_PXYZE.E());
  double a = (P_long * a1)/(2 * a2);
  double r = (((Pvis_PXYZE.E() * Pvis_PXYZE.E()) * a1 * a1)/(4 * (a2* a2))) + (((Pvis_PXYZE.E() * P_trans) * (Pvis_PXYZE.E() * P_trans))/a2);
  if (isnan(TMath::Sqrt(r))){
    r = 0.0;
  }
double b0_P_rec_plus = P_long - a  + TMath::Sqrt(r);
  double b0_P_rec_minus = P_long - a  - TMath::Sqrt(r);
  double b0_P_rec;
  if (sign == 1){
    b0_P_rec = b0_P_rec_plus;
  }
  else{
    b0_P_rec = b0_P_rec_minus;
  }
  double b0_PX_rec = b0_cosangle_X * b0_P_rec;
  double b0_PY_rec = b0_cosangle_Y * b0_P_rec;
  double b0_PZ_rec = b0_cosangle_Z * b0_P_rec;
  double b0_PE_rec = TMath::Sqrt((b0_M_val * b0_M_val) + (b0_P_rec * b0_P_rec));
  TLorentzVector PXYZE_rec(b0_PX_rec, b0_PY_rec, b0_PZ_rec, b0_PE_rec);

  return (PXYZE_rec);
}

double get_q2(TLorentzVector Pcomb)
{
  double q2_p;
  q2_p = Pcomb.E()*Pcomb.E() - Pcomb.X()*Pcomb.X() - Pcomb.Y()*Pcomb.Y() - Pcomb.Z()*Pcomb.Z();

  return (q2_p/1000000);
}

double get_w(double q2)
{
  double w;
  double b0_M_val = 5279.610;
  double dst_M_val = 2010.27;
  double mass_Dst_GeV = dst_M_val / 1000.0;
  double mass_B0_GeV = b0_M_val / 1000.0;
  w = (mass_B0_GeV*mass_B0_GeV + mass_Dst_GeV*mass_Dst_GeV - q2) / (2*mass_B0_GeV*mass_Dst_GeV);
  return w;
}

double get_mmiss2(TLorentzVector Pcomb2, TLorentzVector mu_PXYZE, TLorentzVector dst_PXYZE, TLorentzVector b0_PXYZE)
{
  double mmiss2_p;
  mmiss2_p = Pcomb2.E()*Pcomb2.E() - (Pcomb2.X()*Pcomb2.X()) - (Pcomb2.Y()*Pcomb2.Y()) - (Pcomb2.Z()*Pcomb2.Z());
  return (mmiss2_p/1000000);
}

#include <random>
int get_rand()
{    
    int sign_quad = rand() % 2; //randint(0,1)
    return sign_quad;
}                                                                                                                                                                                              
float calc_eta(float mu_P, float mu_PZ){                                                                                                                                                     
    float eta = 0.5*TMath::Log((mu_P + mu_PZ)/(mu_P - mu_PZ));                                                                                                                               
return eta;                                                                                                                                                                                   
}

TVector3 unit_Z(0,0,1);

double b0_M = 5279.61;
double mu_M = 105.6583745;

TH3D * hist_pi_P_smearing;
TH1D * hpiP_mu_unrolled_cumulative;
TH3D * hist_K_P_smearing;
TH1D * hKP_mu_unrolled_cumulative;
 
void get_smearing_histograms(std::string Pi_smearing_file, std::string K_smearing_file){

TFile *Pi_file = TFile::Open((Pi_smearing_file).c_str(), "READ");
TFile *K_file = TFile::Open((K_smearing_file).c_str(), "READ");
 
hist_pi_P_smearing = (TH3D*)Pi_file->Get("hpiP_mu_diff");
hpiP_mu_unrolled_cumulative = (TH1D*)Pi_file->Get("hpiP_mu_unrolled_cumulative");
 
hist_K_P_smearing = (TH3D*)K_file->Get("hKP_mu_diff");
hKP_mu_unrolled_cumulative = (TH1D*)K_file->Get("hKP_mu_unrolled_cumulative");
return;
}

TRandom3 myrnd;

std::array<double,3> compute_mu_P_corr(double mu_PX, double mu_PY, double mu_PZ, float w_pi, float w_K, float w_p, float w_e, float w_g) {

    std::array<double,3> new_mu_P;
    myrnd.SetSeed(0);
 
    float is_smeared = myrnd.Uniform();
    float fraction_of_smeared = 1.0;
 
    float new_momentum_pion_x;
    float new_momentum_Kaon_x;
    float new_momentum_pion_y;
    float new_momentum_Kaon_y;
    float new_momentum_pion_z;
    float new_momentum_Kaon_z;

    double x;
    double y;
    double z;
 
    double x_Kaon;
    double y_Kaon;
    double z_Kaon;
 
    hist_pi_P_smearing->GetRandom3(x,y,z);

    new_momentum_pion_x = mu_PX + x*mu_PX;
    new_momentum_pion_y = mu_PY + y*mu_PY;
    new_momentum_pion_z = mu_PZ + z*mu_PZ;
 
    hist_K_P_smearing->GetRandom3(x_Kaon,y_Kaon,z_Kaon);
 
    new_momentum_Kaon_x = mu_PX + x_Kaon*mu_PX;
    new_momentum_Kaon_y = mu_PY + y_Kaon*mu_PY;
    new_momentum_Kaon_z = mu_PZ + z_Kaon*mu_PZ;
 
    std::array<double,3> new_mu_P_pion{new_momentum_pion_x, new_momentum_pion_y, new_momentum_pion_z};
    std::array<double,3> new_mu_P_Kaon{new_momentum_Kaon_x, new_momentum_Kaon_y, new_momentum_Kaon_z};
    std::array<double,3> mu_P{mu_PX, mu_PY, mu_PZ};

    if((w_pi>w_K)&&(w_pi>w_g)&&(w_pi>w_p)&&(w_pi>w_e)){
        if(is_smeared<fraction_of_smeared){
            new_mu_P = new_mu_P_pion;
        }
        else{
            new_mu_P = mu_P;
        }
    }
    else if((w_K>w_pi)&&(w_K>w_g)&&(w_K>w_p)&&(w_K>w_e)){
        new_mu_P = new_mu_P_Kaon;
    }
    else{
        new_mu_P = mu_P;
    }

    return new_mu_P;
}

TLorentzVector smeared_momentum(std::array<double,3> mu_P_corr){

    TLorentzVector mu_LV;
    mu_LV.SetXYZM(mu_P_corr[0],mu_P_corr[1],mu_P_corr[2],mu_M);
    return mu_LV;
}

'''

ROOT.gInterpreter.Declare(auxiliar_functions_code)

tree_branches = [i.GetName() for i in tin.GetListOfBranches()]

branches_to_forget = []
for br in branches_to_forget:
    tin.SetBranchStatus(br,0)
    print("branch OFF ", br)

df = ROOT.RDataFrame(tin)
df_base = df.Filter("pass_selection > 0") # base selections like Mass ,IP, BDT, P, PT
df22 = df_base.Define("b0_PXYZE", "TLorentzVector(b0_PX, b0_PY, b0_PZ, b0_PE)") \
         .Define("mu_PXYZE", "TLorentzVector(mu_PX, mu_PY, mu_PZ, mu_PE)") \
         .Define("dst_PXYZE", "TLorentzVector(dst_PX, dst_PY, dst_PZ, dst_PE)") \
         .Define("d0_PXYZE", "TLorentzVector(d0_PX, d0_PY, d0_PZ, d0_PE)") \
         .Define("pi_PXYZE", "TLorentzVector(spi_PX, spi_PY, spi_PZ, spi_PE)") \
         .Define("b0_ENDVERTEX_XYZ", "TVector3(b0_ENDVERTEX_X, b0_ENDVERTEX_Y, b0_ENDVERTEX_Z)") \
         .Define("b0_OWNPV_XYZ", "TVector3(b0_OWNPV_X, b0_OWNPV_Y, b0_OWNPV_Z)") \
         .Define("b0_DIRECTION", "TVector3(b0_ENDVERTEX_XYZ.X() - b0_OWNPV_XYZ.X(), b0_ENDVERTEX_XYZ.Y() - b0_OWNPV_XYZ.Y(), b0_ENDVERTEX_XYZ.Z() - b0_OWNPV_XYZ.Z())") \
         .Define("Pvis_PXYZE", "TLorentzVector(mu_PXYZE.X() + dst_PXYZE.X(), mu_PXYZE.Y() + dst_PXYZE.Y(), mu_PXYZE.Z() + dst_PXYZE.Z(), mu_PXYZE.E() + dst_PXYZE.E())") \
         .Define("b0_PXYZE_rec_Brestframe", "get_Bmom_rec_Brestframe(b0_DIRECTION, dst_PXYZE, mu_PXYZE, b0_M, b0_PZ)") \
         .Define("nu_PXYZE_rec_Brestframe", "TLorentzVector(b0_PXYZE_rec_Brestframe.X() - dst_PXYZE.X() - mu_PXYZE.X(), b0_PXYZE_rec_Brestframe.Y() - dst_PXYZE.Y() - mu_PXYZE.Y(), b0_PXYZE_rec_Brestframe.Z() - dst_PXYZE.Z() - mu_PXYZE.Z(), b0_PXYZE_rec_Brestframe.E() - dst_PXYZE.E() - mu_PXYZE.E())") \
         .Define("W_PXYZE_rec_Brestframe", "TLorentzVector(b0_PXYZE_rec_Brestframe.X() - dst_PXYZE.X(), b0_PXYZE_rec_Brestframe.Y() - dst_PXYZE.Y(), b0_PXYZE_rec_Brestframe.Z() - dst_PXYZE.Z(), b0_PXYZE_rec_Brestframe.E() - dst_PXYZE.E())") \
         .Define("sign_quad", "get_rand()") \
         .Define("b0_PXYZE_rec_quadratic", "get_Bmom_rec_quadratic(b0_DIRECTION, dst_PXYZE, mu_PXYZE, Pvis_PXYZE, sign_quad)") \
         .Define("b0_PXYZE_rec_quadratic_plus", "get_Bmom_rec_quadratic(b0_DIRECTION, dst_PXYZE, mu_PXYZE, Pvis_PXYZE, 1)") \
         .Define("b0_PXYZE_rec_quadratic_minus", "get_Bmom_rec_quadratic(b0_DIRECTION, dst_PXYZE, mu_PXYZE, Pvis_PXYZE, 0)") \
         .Define("nu_PXYZE_rec_quadratic", "TLorentzVector(b0_PXYZE_rec_quadratic.X() - dst_PXYZE.X() - mu_PXYZE.X(), b0_PXYZE_rec_quadratic.Y() - dst_PXYZE.Y() - mu_PXYZE.Y(), b0_PXYZE_rec_quadratic.Z() - dst_PXYZE.Z() - mu_PXYZE.Z(), b0_PXYZE_rec_quadratic.E() - dst_PXYZE.E() - mu_PXYZE.E())") \
         .Define("W_PXYZE_rec_quadratic", "TLorentzVector(b0_PXYZE_rec_quadratic.X() - dst_PXYZE.X(), b0_PXYZE_rec_quadratic.Y() - dst_PXYZE.Y(), b0_PXYZE_rec_quadratic.Z() - dst_PXYZE.Z(), b0_PXYZE_rec_quadratic.E() - dst_PXYZE.E())") \
         .Define("Pcomb0_Brestframe", "TLorentzVector(b0_PXYZE_rec_Brestframe.X() - dst_PXYZE.X(), b0_PXYZE_rec_Brestframe.Y() - dst_PXYZE.Y(), b0_PXYZE_rec_Brestframe.Z() - dst_PXYZE.Z(), b0_PXYZE_rec_Brestframe.E() - dst_PXYZE.E())") \
         .Define("Pcomb0_quadratic", "TLorentzVector(b0_PXYZE_rec_quadratic.X() - dst_PXYZE.X(), b0_PXYZE_rec_quadratic.Y() - dst_PXYZE.Y(), b0_PXYZE_rec_quadratic.Z() - dst_PXYZE.Z(), b0_PXYZE_rec_quadratic.E() - dst_PXYZE.E())") \
         .Define("Pcomb0_quadratic_plus", "TLorentzVector(b0_PXYZE_rec_quadratic_plus.X() - dst_PXYZE.X(), b0_PXYZE_rec_quadratic_plus.Y() - dst_PXYZE.Y(), b0_PXYZE_rec_quadratic_plus.Z() - dst_PXYZE.Z(), b0_PXYZE_rec_quadratic_plus.E() - dst_PXYZE.E())") \
         .Define("Pcomb0_quadratic_minus", "TLorentzVector(b0_PXYZE_rec_quadratic_minus.X() - dst_PXYZE.X(), b0_PXYZE_rec_quadratic_minus.Y() - dst_PXYZE.Y(), b0_PXYZE_rec_quadratic_minus.Z() - dst_PXYZE.Z(), b0_PXYZE_rec_quadratic_minus.E() - dst_PXYZE.E())") \
         .Define("Pcomb2_Brestframe", "TLorentzVector(b0_PXYZE_rec_Brestframe.X() - dst_PXYZE.X() - mu_PXYZE.X(), b0_PXYZE_rec_Brestframe.Y() - dst_PXYZE.Y() - mu_PXYZE.Y(), b0_PXYZE_rec_Brestframe.Z()- dst_PXYZE.Z() - mu_PXYZE.Z(), b0_PXYZE_rec_Brestframe.E() - dst_PXYZE.E() - mu_PXYZE.E())") \
         .Define("Pcomb2_quadratic", "TLorentzVector(b0_PXYZE_rec_quadratic.X() - dst_PXYZE.X() - mu_PXYZE.X(), b0_PXYZE_rec_quadratic.Y() - dst_PXYZE.Y() - mu_PXYZE.Y(), b0_PXYZE_rec_quadratic.Z()- dst_PXYZE.Z() - mu_PXYZE.Z(), b0_PXYZE_rec_quadratic.E() - dst_PXYZE.E() - mu_PXYZE.E())") \

df23 = df22.Define("nu_PE","nu_PXYZE_rec_quadratic.E()")
df24 = df23.Define("nu_PX","nu_PXYZE_rec_quadratic.X()")
df25 = df24.Define("nu_PY","nu_PXYZE_rec_quadratic.Y()")
df26 = df25.Define("nu_PZ","nu_PXYZE_rec_quadratic.Z()")
df27 = df26.Define("costhetalrec_brf", "get_costhetal(b0_PXYZE_rec_Brestframe,mu_PXYZE,W_PXYZE_rec_Brestframe)")
df28 = df27.Define("costhetadrec_brf", "get_costhetad(b0_PXYZE_rec_Brestframe,d0_PXYZE,dst_PXYZE,W_PXYZE_rec_Brestframe)")
df29 = df28.Define("chirec_brf", "get_chi(nu_PXYZE_rec_Brestframe,mu_PXYZE,pi_PXYZE,d0_PXYZE,b0_PXYZE_rec_Brestframe,W_PXYZE_rec_Brestframe,dst_PXYZE)")
df30 = df29.Define("costhetalrec_quad", "get_costhetal(b0_PXYZE_rec_quadratic,mu_PXYZE,W_PXYZE_rec_quadratic)")
df31 = df30.Define("costhetadrec_quad", "get_costhetad(b0_PXYZE_rec_quadratic,d0_PXYZE,dst_PXYZE,W_PXYZE_rec_quadratic)")
df32 = df31.Define("chirec_quad", "get_chi(nu_PXYZE_rec_quadratic,mu_PXYZE,pi_PXYZE,d0_PXYZE,b0_PXYZE_rec_quadratic,W_PXYZE_rec_quadratic,dst_PXYZE)")
df33 = df32.Define("q2rec_brf", "get_q2(Pcomb0_Brestframe)")
df34 = df33.Define("q2rec_quad", "get_q2(Pcomb0_quadratic)")
df35 = df34.Define("mmiss2rec_brf", "get_mmiss2(Pcomb2_Brestframe,mu_PXYZE,dst_PXYZE,b0_PXYZE_rec_Brestframe)")
df36 = df35.Define("mmiss2rec_quad", "get_mmiss2(Pcomb2_quadratic,mu_PXYZE,dst_PXYZE,b0_PXYZE_rec_quadratic)")
df37 = df36.Define("deltamass", "dst_M - d0_M")
df38 = df37.Define("b0_Prec_quad", "b0_PXYZE_rec_quadratic.P()")
df39 = df38.Define("b0_Prec_brf", "b0_PXYZE_rec_Brestframe.P()")
df40 = df39.Define("w_rec_quad", "get_w(q2rec_quad)")
df41 = df40.Define("q2rec_plus", "get_q2(Pcomb0_quadratic_plus)")
df42 = df41.Define("q2rec_minus", "get_q2(Pcomb0_quadratic_minus)")
df43 = df42.Define("b0_inv_mass" , "b0_PXYZE.M()")
df44 = df43.Define("dst_inv_mass" , "dst_PXYZE.M()")
# Add weight branches
df45 = df44.Define("pid_eff", "PIDCalibEff")

         #.Define("b0_ENDVERTEX_XYZ_TRUE", "TVector3(b0_TRUEENDVERTEX_X, b0_TRUEENDVERTEX_Y, b0_TRUEENDVERTEX_Z)") \
         #.Define("b0_DIRECTION_TRUE", "TVector3(b0_ENDVERTEX_XYZ_TRUE.X() - b0_OWNPV_XYZ_TRUE.X(), b0_ENDVERTEX_XYZ_TRUE.Y() - b0_OWNPV_XYZ_TRUE.Y(), b0_ENDVERTEX_XYZ_TRUE.Z() - b0_OWNPV_XYZ_TRUE.Z())")\
         #.Define("b0_OWNPV_XYZ_TRUE", "TVector3(b0_TRUEORIGINVERTEX_X, b0_TRUEORIGINVERTEX_Y, b0_TRUEORIGINVERTEX_Z)") \

df_truth1 = df45.Define("b0_PXYZE_TRUE", "TLorentzVector(b0_TRUEP_X, b0_TRUEP_Y, b0_TRUEP_Z, b0_TRUEP_E)") \
         .Define("mu_PXYZE_TRUE", "TLorentzVector(mu_TRUEP_X, mu_TRUEP_Y, mu_TRUEP_Z, mu_TRUEP_E)") \
         .Define("dst_PXYZE_TRUE", "TLorentzVector(dst_TRUEP_X, dst_TRUEP_Y, dst_TRUEP_Z, dst_TRUEP_E)") \
         .Define("d0_PXYZE_TRUE", "TLorentzVector(d0_TRUEP_X, d0_TRUEP_Y, d0_TRUEP_Z, d0_TRUEP_E)") \
         .Define("pi_PXYZE_TRUE", "TLorentzVector(spi_TRUEP_X, spi_TRUEP_Y, spi_TRUEP_Z, spi_TRUEP_E)") \
         .Define("nu_PXYZE_TRUE", "TLorentzVector(b0_PXYZE_TRUE.X() - dst_PXYZE_TRUE.X() - mu_PXYZE_TRUE.X(), b0_PXYZE_TRUE.Y() - dst_PXYZE_TRUE.Y() - mu_PXYZE_TRUE.Y(), b0_PXYZE_TRUE.Z() - dst_PXYZE_TRUE.Z() - mu_PXYZE_TRUE.Z(), b0_PXYZE_TRUE.E() - dst_PXYZE_TRUE.E() - mu_PXYZE_TRUE.E())") \
         .Define("W_PXYZE_TRUE", "TLorentzVector(b0_PXYZE_TRUE.X() - dst_PXYZE_TRUE.X(), b0_PXYZE_TRUE.Y() - dst_PXYZE_TRUE.Y(), b0_PXYZE_TRUE.Z() - dst_PXYZE_TRUE.Z(), b0_PXYZE_TRUE.E() - dst_PXYZE_TRUE.E())") \
         .Define("Pcomb0_TRUE", "TLorentzVector(b0_PXYZE_TRUE.X() - dst_PXYZE_TRUE.X(), b0_PXYZE_TRUE.Y() - dst_PXYZE_TRUE.Y(), b0_PXYZE_TRUE.Z() - dst_PXYZE_TRUE.Z(), b0_PXYZE_TRUE.E() - dst_PXYZE_TRUE.E())") \
         .Define("Pcomb2_TRUE", "TLorentzVector(b0_PXYZE_TRUE.X() - dst_PXYZE_TRUE.X() - mu_PXYZE_TRUE.X(), b0_PXYZE_TRUE.Y() - dst_PXYZE_TRUE.Y() - mu_PXYZE_TRUE.Y(), b0_PXYZE_TRUE.Z() - dst_PXYZE_TRUE.Z() - mu_PXYZE_TRUE.Z(), b0_PXYZE_TRUE.E() - dst_PXYZE_TRUE.E() - mu_PXYZE_TRUE.E())") \
 
df_truth2 = df_truth1.Define("costhetal_TRUE", "get_costhetal(b0_PXYZE_TRUE,mu_PXYZE_TRUE,W_PXYZE_TRUE)")
df_truth3 = df_truth2.Define("costhetad_TRUE", "get_costhetad(b0_PXYZE_TRUE,d0_PXYZE_TRUE,dst_PXYZE_TRUE,W_PXYZE_TRUE)")
df_truth4 = df_truth3.Define("chi_TRUE", "get_chi(nu_PXYZE_TRUE,mu_PXYZE_TRUE,pi_PXYZE_TRUE,d0_PXYZE_TRUE,b0_PXYZE_TRUE,W_PXYZE_TRUE,dst_PXYZE_TRUE)")
df_truth5 = df_truth4.Define("q2_TRUE", "get_q2(Pcomb0_TRUE)")
df_truth6 = df_truth5.Define("mmiss2_TRUE", "get_mmiss2(Pcomb2_TRUE,mu_PXYZE_TRUE,dst_PXYZE_TRUE,b0_PXYZE_TRUE)")
df_truth9 = df_truth6.Define("b0_P_TRUE", "b0_PXYZE_TRUE.P()")
df_truth10 = df_truth9.Define("nu_TRUEP_E","nu_PXYZE_TRUE.E()")
df_truth11 = df_truth10.Define("nu_TRUEP_X","nu_PXYZE_TRUE.X()")
df_truth12 = df_truth11.Define("nu_TRUEP_Y","nu_PXYZE_TRUE.Y()")
df_truth13 = df_truth12.Define("nu_TRUEP_Z","nu_PXYZE_TRUE.Z()")
df_truth14 = df_truth13.Define("w_TRUE", "get_w(q2_TRUE)")
df_truth15 = df_truth14.Define("dst_M_TRUE", "dst_PXYZE_TRUE.M()")
df_truth16 = df_truth15.Define("b0_M_TRUE", "b0_PXYZE_TRUE.M()")

column_name_vector = ROOT.std.vector('string')()
column_name_vector.push_back("q2rec_quad")
column_name_vector.push_back("costhetalrec_quad")
column_name_vector.push_back("costhetadrec_quad")
column_name_vector.push_back("chirec_quad")
column_name_vector.push_back("mmiss2rec_quad")
column_name_vector.push_back("w_rec_quad")
column_name_vector.push_back("q2rec_plus")
column_name_vector.push_back("q2rec_minus")
column_name_vector.push_back("b0_Prec_quad")
column_name_vector.push_back("b0_Prec_brf")
column_name_vector.push_back("nu_PE")
column_name_vector.push_back("nu_PX")
column_name_vector.push_back("nu_PY")
column_name_vector.push_back("nu_PZ")
column_name_vector.push_back("b0_inv_mass")
column_name_vector.push_back("dst_inv_mass")
column_name_vector.push_back("pid_eff")

column_name_vector.push_back("q2rec_brf")
column_name_vector.push_back("costhetalrec_brf")
column_name_vector.push_back("costhetadrec_brf")
column_name_vector.push_back("chirec_brf")
column_name_vector.push_back("mmiss2rec_brf")

column_name_vector.push_back("q2_TRUE")
column_name_vector.push_back("costhetal_TRUE")
column_name_vector.push_back("costhetad_TRUE")
column_name_vector.push_back("chi_TRUE")
column_name_vector.push_back("mmiss2_TRUE")
column_name_vector.push_back("nu_TRUEP_E")
column_name_vector.push_back("nu_TRUEP_X")
column_name_vector.push_back("nu_TRUEP_Y")
column_name_vector.push_back("nu_TRUEP_Z")
column_name_vector.push_back("w_TRUE")
column_name_vector.push_back("b0_M_TRUE")
column_name_vector.push_back("dst_M_TRUE")

for br in tree_branches:
    if not (br in branches_to_forget):
        column_name_vector.push_back(br)

df_truth16.Snapshot("DecayTree", outfilename, column_name_vector)
print("Finished writing to file ", outfilename)