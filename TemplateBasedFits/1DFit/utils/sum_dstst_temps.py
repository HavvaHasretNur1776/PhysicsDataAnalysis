import numpy as np
import os

# rel_BFs
# dec_fracs
# gen_yields
rel_BFs = {
 
    'b_dstdx' : {
        "b+_dstdx": [0.5, 0],
        "b0_dstdx": [0.5, 0],
    },
    'b_dststmu' : {
        "b0_dststmu_415": [0.5 * 6.6 * 1e-4, 0.5 * 1.1 * 1e-4],#0.5 factors for D**- -> D*- isospin symmetry, keep original measurement errors
        "b0_dststmu_10413": [0.5 * 2.85 * 1e-3, 0.5 * 0.25*1e-3],
        "b0_dststmu_20413": [0.5 * 2.5 *1e-3, 0.5 * 0.5 * 1e-3],
        
        "b+_dststmu_425": [1.06*1e-3, 0.28*1e-3],
        "b+_dststmu_10423": [2.84 * 1e-3, 0.17 * 1e-3],
        "b+_dststmu_20423": [1.7 * 1e-3, 0.6 * 1e-3],
    },
    'bs_dststmu' : {
        "bs_dststmu_10433": [2 * 2.7 * 1e-3, 2 * 0.7 * 1e-3],
        "bs_dststmu_435": [4.3 * 1e-5, 2.24 * 1e-5],

    },
    'b_dststmu_heavy' : {
        "b+_dststmu_heavy": [0.5 * 6.6 * 1e-4, 0.5 * 1.1 * 1e-4],
        "b0_dststmu_heavy": [1.06 * 1e-3, 0.18 * 1e-3],
    },
}



dec_fracs = {
    'b_dstdx' : {
        "b+_dstdx": 1,
        "b0_dstdx": 1,
    },
    'b_dststmu' : {
        "b0_dststmu_415": 0.1188 * (0.19 + 0.0008),
        "b0_dststmu_10413": 0.2727 * 0.1829,
        "b0_dststmu_20413": 0.3069 * 0.4547,
        
        "b+_dststmu_425": 0.1471 * (0.4750 + 0.0014),
        "b+_dststmu_10423": 0.2341 * 0.3821,
        "b+_dststmu_20423": 0.2647 * 0.3821,
    },
    'bs_dststmu' : {
        "bs_dststmu_10433": 0.0070 * 0.5,
        "bs_dststmu_435": 0.0070 * 0.5,

    },
    'b_dststmu_heavy' : {
        "b+_dststmu_heavy": 0.05 * 0.08,
        "b0_dststmu_heavy": 0.05 * 0.08,
    },
}


gen_yields = {#2016 (down+up) + 2017 (down+up)+ 2018 (down+up)
    'b_dstdx' : {#eventID=11894610, 12895400 
        "b+_dstdx":(1600224718 + 1595754180) + (1310581403 +1315729955)+(1257034067+1315729955),
        "b0_dstdx": (1256782128 + 1271630284) + (1429806082 +1378395336)+(1021569152+1378395336)
    },
    'b_dststmu' : {#11874430, 12873450
        "b0_dststmu_415": (1199167539 + 1236752126)+(1138905547+1131098996)+(1341232837+1131098996),
        "b0_dststmu_10413": (1199167539 + 1236752126)+(1138905547+1131098996)+(1341232837+1131098996),
        "b0_dststmu_20413":(1199167539 + 1236752126)+(1138905547+1131098996)+(1341232837+1131098996),
        
        "b+_dststmu_425": (1227691430 + 1175883049)+(1140146548+1145355006)+(1301178439+1145355006),
        "b+_dststmu_10423": (1227691430 + 1175883049)+(1140146548+1145355006)+(1301178439+1145355006),
        "b+_dststmu_20423": (1227691430 + 1175883049)+(1140146548+1145355006)+(1301178439+1145355006)
    },
    'bs_dststmu' : {#13674000
        "bs_dststmu_10433": (1492068605 + 1453243128) +(1359828427+1368169830)+(1186885466+1368169830),
        "bs_dststmu_435": (1492068605 + 1453243128)+(1359828427+1368169830)+(1186885466+1368169830)

    },
    'b_dststmu_heavy' : {#11676012, 12675402
        "b+_dststmu_heavy": (1338123297 + 1448221903) + (1218633361+1235141002)+(1133194996+1235141002),
        "b0_dststmu_heavy": (1188664181 + 1198834296)+ (1256688197+1200763722)+(1050410052+1200763722)
    },
}




input_dir = "./"     # your .dat files
output_dir = "grouped_templates"
os.makedirs(output_dir, exist_ok=True)

for group in rel_BFs.keys():

    components = rel_BFs[group].keys()

    # find available wbins from first component
    example_files = [f for f in os.listdir(input_dir)
                     if list(components)[0] in f]

    for file in example_files:

        wbin = file.split("_wbin")[-1].replace(".dat","")

        combined = None

        values_list = []
        physics_weights = []

        for comp in components:

            filepath = os.path.join(input_dir,
                                    f"{comp}_wbin{wbin}.dat")

            if not os.path.exists(filepath):
                continue

            values = np.loadtxt(filepath)

            BF = rel_BFs[group][comp][0]
            dec = dec_fracs[group][comp]
            Ngen = gen_yields[group][comp]

            eff = np.sum(values) / Ngen
            weight = BF * dec * eff

            values_list.append(values)
            physics_weights.append(weight)

        # normalize physics weights
        physics_weights = np.array(physics_weights)
        physics_weights /= np.sum(physics_weights)

        # build grouped histogram
        grouped_hist = sum(w * v for w, v in zip(physics_weights, values_list))

        # restore total yield
        total_naive = sum(np.sum(v) for v in values_list)
        grouped_hist *= total_naive / np.sum(grouped_hist)


        outpath = os.path.join(output_dir, f"{group}_wbin{wbin}.dat")

        np.savetxt(outpath, grouped_hist)

        print(f"Written {outpath}")