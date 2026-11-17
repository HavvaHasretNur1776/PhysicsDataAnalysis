import numpy as np


template_path = '/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/templates_1D/'
template_path_5D = '/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/templates_5D/'
template_path_control = '/eos/lhcb/wg/semileptonic/B2Dstmunu_angular_Run2/templates_control/'
template_paths = {'1D': template_path,
                  '5D': template_path_5D}
tuple_path =  '/data/lhcb/users/martel/shape_measurement/MC/'
data_tuple_path =  '/data/lhcb/users/martel/shape_measurement/data/'


jsons_dir = '/home/martel/LHCb/b02dstmunu_q2_shape/jsons/'

m_corr_bins = 10
n_w_bins = 5
component_names = ["nDstMu", "nDstTau", "nDstDX", "nbDstst", "nbsDstst", "nDstst_higher"]

component_labels = {
    "nDstMu": r"$B^0 \to D^{*+} \mu^- \nu_\mu$",
    "nDstTau": r"$B^0 \to D^{*+} \tau^- \nu_\tau$",
    "nDstDX": r"$B \to D^{*}D(X_c \to \mu \nu X)$",
    "nbDstst": r"$B \to D^{**} \mu^- \nu_\mu$",
    "nbsDstst": r"$B_s \to D^{**} \mu^- \nu_\mu$",
    "nDstst_higher": r"$B \to D^{**}(higher) \mu^- \nu_\mu$",
    'nWSMu' : 'WS Mu',
   # 'nWSPi' : 'WS Pi',
    'data' : 'data'
}

component_colors = {
 'nbsDstst': '#de6dda',
 'nDstTau': '#B8FFCC',
 'nDstDX': (0.4707222689300699, 0.4984208187438168, 0.21398565629427133),
 'nbDstst': '#690065',
'nDstst_higher': '#ba8db9',
 'nDstMu': '#FF0051',
 'nWSMu' : 'gold',
 #'nWSPi' : '#cf5c04',
 'data' : 'black'
 }

bin_map_5D ={}
bin_labels_5D = {}

component_templates = {
 'nbsDstst': 'bs_Dstst_bin',
 'nDstTau': 'b0_DstTau_bin',
 'nDstDX': 'DstDx_bin',
 'nbDstst': 'b_Dstst_bin',
'nDstst_higher': 'Dstst_higher_bin',
 'nWSMu':  'data_WSMu_bin',
 #'nWSPi':  'data_WSPi_bin',
 'nDstMu':  'signal_bin',
 }

mc_templates = {
 'nbsDstst': 'bs_Dstst_bin',
 'nDstTau': 'b0_DstTau_bin',
 'nDstDX': 'DstDx_bin',
 'nbDstst': 'b_Dstst_bin',
'nDstst_higher': 'Dstst_higher_bin',
 'nDstMu':  'signal_bin',
 }

data_path_up = "data_up_bin"
data_path_down =  "data_down_bin"
data_path = "data_bin"

fit_init_params = {
        'scale_nDstMu' : [0.75],
        'scale_nbDstst' : [0.05],
        'scale_nbsDstst' : [0.03],
        'scale_nDstst_higher' : [0.07],
        'scale_nDstTau' : [0.03],
        'scale_nDstDX' : [0.002],
        'scale_nWSMu' : [0.025],
     #   'scale_nWSPi' : [0.025]
    }

fit_fixed_fracs = {
        #'nDstMu' : [0.88],
        'nbDstst' : 0.04,#0.1,
        'nbsDstst' : 0.02,#0.07,
        'nDstst_higher' : 0.01,#0.05,
        'nDstTau' : 0.02,#0.05,
       'nDstDX' : 0.01,#0.04,
        'nWSMu' : 0.05,
       # 'nWSPi' : 0.007
    }

mmiss2_bin_edges = {}
mmiss2_bin_edges['1D_5'] =                    { 
    0: np.linspace(-5.0, 6.0, m_corr_bins + 1),  
    1: np.linspace(-3.0, 5.5, m_corr_bins + 1), 
    2: np.linspace(-3.0, 5.0, m_corr_bins + 1),
    3: np.linspace(-2.8, 5.0, m_corr_bins + 1),
    4: np.linspace(-2.0, 4.1, m_corr_bins + 1)
                    } 


mmiss2_bin_edges['1D_10'] =                    { 
    0: np.linspace(-4.0, 5.5, m_corr_bins + 1),  
    1: np.linspace(-5.0, 6.0, m_corr_bins + 1), 
    2: np.linspace(-3.0, 5.5, m_corr_bins + 1),
    3: np.linspace(-3, 5.0, m_corr_bins + 1),
    4: np.linspace(-3.0, 5.0, m_corr_bins + 1),
    5: np.linspace(-3.0, 5.0, m_corr_bins + 1),  
    6: np.linspace(-2.8, 5.0, m_corr_bins + 1), 
    7: np.linspace(-2.0, 5.0, m_corr_bins + 1),
    8: np.linspace(-2.0, 4.1, m_corr_bins + 1),
    9: np.linspace(-2.0, 3.0, m_corr_bins + 1),
                    } 



mmiss2_bin_edges['control'] = { 
    0: np.linspace(-2.0, 8.0, m_corr_bins + 1),  
    1: np.linspace(-2.0, 8.0, m_corr_bins + 1), 
    2: np.linspace(-2.0, 8.0, m_corr_bins + 1),
    3: np.linspace(-2.0, 8.0, m_corr_bins + 1),
    4: np.linspace(-2.0, 8.0, m_corr_bins + 1)
                    } 

mmiss2_bin_edges_full = np.linspace(-5.0, 6.0, m_corr_bins + 1)

w_edges = {
    5 :  [1.  , 1.1 , 1.2,  1.3 , 1.4 ,  1.51 ],
    10 : [1, 1.050, 1.1,1.150, 1.200, 1.250, 1.300, 1.350, 1.4, 1.450, 1.51]
}
variables_4D = ['w_reco_regreg','ctl_reco_regreg', 'ctd_reco_regreg', 'chi_reco_regreg', 'm2miss_reco_approx']
angular_vars = ['w_reco_regreg','ctl_reco_regreg', 'ctd_reco_regreg', 'chi_reco_regreg']

#w_bins_4D = [1.  , 1.05, 1.1 , 1.15, 1.2 , 1.25, 1.3 , 1.35, 1.4 , 1.45, 1.51 ]
w_bins_4D = [1.  , 1.1 , 1.2,  1.3 , 1.4 ,  1.51 ]

nbins_4D = {'w_reco_regreg': 5,
            'ctl_reco_regreg': 3,
            'ctd_reco_regreg': 3,
            'chi_reco_regreg': 5,
            'm2miss_reco_approx' : 10
            }
bin_no = 0

for bv0 in range(nbins_4D['w_reco_regreg']):
    for bv1 in range(nbins_4D['ctl_reco_regreg']):
        for bv2 in range(nbins_4D['ctd_reco_regreg']):
            for bv3 in range(nbins_4D['chi_reco_regreg']):
                bin_labels_5D[bin_no] = f'$w$ {bv0+1}, $\\cos\\theta_\\ell$ {bv1+1}, $\\cos\\theta_D$ {bv2+1}, $\\chi$ {bv3+1}'
                bin_map_5D[bin_no] = f'{bv0}_{bv1}'
                bin_no+=1
mmiss2_bins_5D = {}
for i in range(nbins_4D['w_reco_regreg']*nbins_4D['ctl_reco_regreg']*nbins_4D['ctd_reco_regreg']*nbins_4D['chi_reco_regreg']):
    if bin_map_5D[i] == '0_0':
        mmiss2_bins_5D[i] = np.linspace(-5, 4, nbins_4D['m2miss_reco_approx'] + 1)
    elif bin_map_5D[i] == '1_0':
        mmiss2_bins_5D[i] = np.linspace(-4, 2.5, nbins_4D['m2miss_reco_approx'] + 1)
    elif bin_map_5D[i] == '1_1':
        mmiss2_bins_5D[i] = np.linspace(-4, 4, nbins_4D['m2miss_reco_approx'] + 1)
    elif bin_map_5D[i] == '2_0':
        mmiss2_bins_5D[i] = np.linspace(-2., 2., nbins_4D['m2miss_reco_approx'] + 1)
    elif bin_map_5D[i] == '2_1':
        mmiss2_bins_5D[i] = np.linspace(-4, 3., nbins_4D['m2miss_reco_approx'] + 1)
    elif bin_map_5D[i] == '3_0':
        mmiss2_bins_5D[i] = np.linspace(-1.5, 1.5, nbins_4D['m2miss_reco_approx'] + 1)
    elif bin_map_5D[i] == '3_1':
        mmiss2_bins_5D[i] = np.linspace(-3., 2.5, nbins_4D['m2miss_reco_approx'] + 1)
    elif bin_map_5D[i] == '4_0':
        mmiss2_bins_5D[i] = np.linspace(-0.5, 1, nbins_4D['m2miss_reco_approx'] + 1)
    elif bin_map_5D[i] == '4_1':
        mmiss2_bins_5D[i] = np.linspace(-1., 2., nbins_4D['m2miss_reco_approx'] + 1)
    else:
        mmiss2_bins_5D[i] = np.linspace(-5, 6, nbins_4D['m2miss_reco_approx'] + 1)

mmiss2_bin_edges['5D_5'] = mmiss2_bins_5D
#mmiss2_bin_edges['5D_5'] = mmiss2_bin_edges['1D_5']
q2_bins = nbins_4D['w_reco_regreg']
ranges_4D = {'w_reco_regreg': (1,1.51),
            'ctl_reco_regreg': (-1,1),
            'ctd_reco_regreg': (-1,1),
            'chi_reco_regreg': (-np.pi,np.pi),
            'm2miss_reco_approx' : (-5.0,6.0)
            }

bin_file = 'binning_hasret_4D_5_3_3_5_10opti_equal_dynsamic_mmiss.json'

gen_yields = {}
gen_yields['2016'] = {
    'nDstDX' : {
        "b+_dstdx":{'up':1595754180,
                     'down':1600224718},
        "b0_dstdx": {'up':1271630284,
                     'down':1256782128}
    },
    'nbDstst' : {
        "b0_dststmu_415": {'up':1236752126,
                     'down':1199167539},
        "b0_dststmu_10413": {'up':1236752126,
                     'down':1199167539}, 
        "b0_dststmu_20413":{'up':1236752126,
                     'down':1199167539}, 
        
        "b+_dststmu_425":{'up':1175883049,
                     'down':1227691430 }, 
        "b+_dststmu_10423":{'up':1175883049,
                     'down':1227691430},  
        "b+_dststmu_20423":{'up':1175883049,
                     'down':1227691430 }, 
    },
    'nbsDstst' : {
        "bs_dststmu_10433":{'up':1453243128,
                     'down':1492068605 }, 
        "bs_dststmu_435":{'up':1453243128,
                     'down':1492068605}, 

    },
    'nDstst_higher' : {
        "b+_dststmu_heavy":{'up':1448221903,
                     'down':1338123297 }, 
        "b0_dststmu_heavy":{'up':1198834296,
                     'down':1188664181 }, 
    },
}

gen_yields['2017'] = {
    'nDstMu' : { 'up' : 1098724106,
                'down' :  1077911703,

    },
    'nDstTau' : { 'up' : 1307755767,
                'down' :  1261331888,

    },

    'nbDstst' : {
        "b0_dststmu_415": {'up':1306807043,
                     'down':1287201281},
        "b0_dststmu_10413": {'up':1306807043,
                     'down':1287201281}, 
        "b0_dststmu_20413":{'up':1306807043,
                     'down':1287201281}, 
        
        "b+_dststmu_425":{'up':1323385459,
                     'down':1277094888}, 
        "b+_dststmu_10423":{'up':1323385459,
                     'down':1277094888},  
        "b+_dststmu_20423":{'up':1323385459,
                     'down':1277094888}, 
    },
    'nbsDstst' : {
        "bs_dststmu_10433":{'up':1322974615,
                     'down':1383608085}, 
        "bs_dststmu_435":{'up':1322974615,
                     'down':1383608085}, 

    },
    'nDstst_higher' : {
        "b+_dststmu_heavy":{'up':1235141002,
                     'down':1218633361}, 
        "b0_dststmu_heavy":{'up':1200763722,
                     'down':1256688197}, 
    },

    'nDstDX' : {
        "b+_dstdx":{'up':1315729955,
                     'down':1310581403},
        "b0_dstdx": {'up':1378395336,
                     'down':1429806082}
    },
}


gen_yields['2018'] = {
    'nDstMu' : { 'up' : 1019733889,
                'down' :  1069054819,

    },
    'nDstTau' : { 'up' : 1094348001,
                'down' :  1055772187,

    },

    'nbDstst' : {
        "b0_dststmu_415": {'up':1139323299,
                     'down':1198013754},
        "b0_dststmu_10413": {'up':1139323299,
                     'down':1198013754}, 
        "b0_dststmu_20413":{'up':1139323299,
                     'down':1198013754}, 
        
        "b+_dststmu_425":{'up':1008915156,
                     'down':1089184362}, 
        "b+_dststmu_10423":{'up':1008915156,
                     'down':1089184362},  
        "b+_dststmu_20423":{'up':1008915156,
                     'down':1089184362}, 
    },
    'nbsDstst' : {
        "bs_dststmu_10433":{'up':1189651648,
                     'down':1219039360}, 
        "bs_dststmu_435":{'up':1189651648,
                     'down':1219039360}, 

    },
    'nDstst_higher' : {
        "b+_dststmu_heavy":{'up':1098571770 ,
                     'down':1133194996}, 
        "b0_dststmu_heavy":{'up':1064433083,
                     'down':1050410052}, 
    },

    'nDstDX' : {
        "b+_dstdx":{'up':1200918897,
                     'down':1257034067},
        "b0_dstdx": {'up':1222293921,
                     'down':1021569152}
    },
}
gen_yields['Run2'] = {
    'nDstMu': {'up': 3235734531, 'down': 3258736817},
 'nDstTau': {'up': 3665280674, 'down': 3535333569},
 'nbDstst': {'b0_dststmu_415': {'up': 3816342157, 'down': 3864389257},
  'b0_dststmu_10413': {'up': 3816342157, 'down': 3864389257},
  'b0_dststmu_20413': {'up': 3816342157, 'down': 3864389257},
  'b+_dststmu_425': {'up': 3678838250, 'down': 3789581347},
  'b+_dststmu_10423': {'up': 3678838250, 'down': 3789581347},
  'b+_dststmu_20423': {'up': 3678838250, 'down': 3789581347}},
 'nbsDstst': {'bs_dststmu_10433': {'up': 3965869391, 'down': 4094716050},
  'bs_dststmu_435': {'up': 3965869391, 'down': 4094716050}},
 'nDstst_higher': {'b+_dststmu_heavy': {'up': 3781934675, 'down': 3689951654},
  'b0_dststmu_heavy': {'up': 3464031101, 'down': 3495762430}},
 'nDstDX': {'b+_dstdx': {'up': 4112403032, 'down': 4167840188},
  'b0_dstdx': {'up': 3872319541, 'down': 3708157362}}
  }

postfit_branch_info = { 
    # B0
    "b0_PT": (100, (0.0, 40000.0), r"$p_{T}(B^{0})\,[\mathrm{MeV}/c]$"),
    "b0_P": (100, (10000.0, 200000.0), r"$p(B^{0})\,[\mathrm{MeV}/c]$"),
    "b0_ETA": (24, (2.0, 6.0), r"$\eta(B^{0})$"),
    "b0_IP_OWNPV": (50, (0.0, 2.0), r"$\mathrm{IP}(B^{0})\,[\mathrm{mm}]$"),
    "b0_FDCHI2_OWNPV": (50, (0.0, 10000.0), r"$\mathrm{FD}\chi^{2}(B^{0})$"),
    "b0_OWNPV_NDOF": (75, (0.0, 300.0), "PV/ndf"),

    # Muon
    "mu_PT": (80, (1000.0, 20000.0), r"$p_{T}(\mu)\,[\mathrm{MeV}/c]$"),
    "mu_P": (100, (10000.0, 100000.0), r"$p(\mu)\,[\mathrm{MeV}/c]$"),
    "mu_ETA": (24, (2.0, 6.0), r"$\eta(\mu)$"),

    # D0
    "d0_PT": (100, (2000.0, 15000.0), r"$p_{T}(D^{0})\,[\mathrm{MeV}/c]$"),
    "d0_P": (80, (10000.0, 150000.0), r"$p(D^{0})\,[\mathrm{MeV}/c]$"),
    "d0_ETA": (24, (2.0, 6.0), r"$\eta(D^{0})$"),
    "d0_IP_OWNPV": (50, (0.0, 2.0), r"$\mathrm{IP}(D^{0})\,[\mathrm{mm}]$"),
    "d0_FDCHI2_OWNPV": (50, (0.0, 100000.0), r"$\mathrm{FD}\chi^{2}(D^{0})$"),

    # D*
    #"dst_FDCHI2_OWNPV": (50, (0.0, 1000.0), r"$\mathrm{FD}\chi^{2}(D^{*})$"),

    # Kaon
    "k_PT": (80, (900.0, 10000.0), r"$p_{T}(K)\,[\mathrm{MeV}/c]$"),
    "k_P": (100, (0.0, 100000.0), r"$p(K)\,[\mathrm{MeV}/c]$"),
    "k_ETA": (24, (2.0, 6.0), r"$\eta(K)$"),
    "k_IP_OWNPV": (50, (0.0, 5.0), r"$\mathrm{IP}(K)\,[\mathrm{mm}]$"),

    # Pion
    "pi_PT": (80, (900.0, 9000.0), r"$p_{T}(\pi)\,[\mathrm{MeV}/c]$"),
    "pi_P": (100, (0.0, 100000.0), r"$p(\pi)\,[\mathrm{MeV}/c]$"),
    "pi_ETA": (24, (2.0, 6.0), r"$\eta(\pi)$"),

    "nTracks": (20, (0.0, 400.0), "nTracks"),
}
