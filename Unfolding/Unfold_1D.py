#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ROOT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import array
import mplhep as hep  # LHCb style
import os
# Apply LHCb style
hep.style.use("LHCb2")

ROOT.gSystem.Load("/eos/lhcb/user/h/havva/q2Measurement-BtoDstMuNu/unfolding/RooUnfold/build/libRooUnfold.so")
ROOT.TH1.AddDirectory(False)

n_iter = 4  # nominal 4, try 3 and 5 for systematics
figdir = "plots"
os.makedirs(figdir, exist_ok=True)

#year = 2018
#input_file_up = f'/eos/lhcb/user/d/dvicoben/bd2dstmunu_run2/Data/MCTracker/{year}/selected/b0_dstmunu_up_selected.root'
#input_file_down = f'/eos/lhcb/user/d/dvicoben/bd2dstmunu_run2/Data/MCTracker/{year}/selected/b0_dstmunu_down_selected.root'

#df_in = ROOT.RDataFrame("DecayTree", str(input_file))
#files = [input_file_up, input_file_down]
#df_in = ROOT.RDataFrame("DecayTree", files)
years = [2016, 2017, 2018]
polarities = ['up', 'down']

# Generate the full list of files
all_files = [
    f'/eos/lhcb/user/d/dvicoben/bd2dstmunu_run2/Data/MCTracker/{y}/selected/b0_dstmunu_{p}_selected.root'
    for y in years for p in polarities
]
df_in = ROOT.RDataFrame("DecayTree", all_files)
df =  pd.DataFrame(df_in.AsNumpy(sorted(['w_reco_regreg', 'w_reco_TRUE'])))
#df = df_.sample(frac = 0.3)these wwere mocking the unfolding 
#df2 = df_.sample(frac = 0.01)this is the data to be unfolded
real_data_files = [
    f'/eos/lhcb/user/d/dvicoben/bd2dstmunu_run2/Data/Collision/{y}/selected/data_RS_{p}_selected.root'
    for y in years for p in polarities
]
df_in2 = ROOT.RDataFrame("DecayTree", real_data_files)
df2 =  pd.DataFrame(df_in2.AsNumpy(sorted(['w_reco_regreg'])))

v_reco = df['w_reco_regreg'].values
v_true = df['w_reco_TRUE'].values
v_reco2 = df2['w_reco_regreg'].values
#v_true2 = df2['w_reco_TRUE'].values

nbins = 10
xmin = 1.0
xmax = 1.5

h_reco = ROOT.TH1F("h_reco", "Reco distribution", nbins, xmin, xmax)
h_true = ROOT.TH1F("h_true", "True distribution", nbins, xmin, xmax)

h_reco2 = ROOT.TH1F("h_reco2", "Reco distribution2", nbins, xmin, xmax)
#h_true2 = ROOT.TH1F("h_true2", "True distribution2", nbins, xmin, xmax)
response = ROOT.RooUnfoldResponse(nbins, xmin, xmax, nbins, xmin, xmax)

for reco, true in zip(v_reco, v_true):
    h_reco.Fill(reco)
    h_true.Fill(true)
    response.Fill(reco, true)

for reco in v_reco2:
    h_reco2.Fill(reco)
#for reco, true in zip(v_reco2, v_true2):
#    h_reco2.Fill(reco)
    #h_true2.Fill(true)



matrix = np.zeros((nbins, nbins))  ######rows = true bins, cols = reco bins
hresponse = response.Hresponse() 
for ix in range(1, nbins+1):
    for iy in range(1, nbins+1):
        matrix[iy-1, ix-1] = hresponse.GetBinContent(ix, iy)

plt.figure(figsize=(8, 6))
im = plt.imshow(matrix, origin='lower', aspect='auto', interpolation='none',  cmap='RdPu')

plt.colorbar(im, label='Counts')
plt.xlabel('Reconstructed bins')
plt.ylabel('True bins')
plt.savefig(f'{figdir}/response_matrix_w_{nbins}_{n_iter}.png', dpi=300)
plt.close()

print(response)
print(h_reco2)
unfold = ROOT.RooUnfoldBayes(response, h_reco2, n_iter)
h_unfolded = unfold.Hunfold()
print('bin err')
print(h_unfolded.GetBinError(1)) 


nbins = h_unfolded.GetNbinsX()
bin_centers = np.array([h_unfolded.GetBinCenter(i) for i in range(1, nbins+1)])
contents = np.array([h_unfolded.GetBinContent(i) for i in range(1, nbins+1)])
errors = np.array([h_unfolded.GetBinError(i) for i in range(1, nbins+1)])
#true_p = np.array([h_true2.GetBinContent(i) for i in range(1, nbins+1)])
reco_p = np.array([h_reco2.GetBinContent(i) for i in range(1, nbins+1)])

bin_width = bin_centers[1] - bin_centers[0]
fig, ax1 = plt.subplots(figsize=(10, 6))
#pulls = (  contents- true_p)/np.sqrt( true_p + errors**2)
#fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True ,gridspec_kw={'height_ratios': [3, 0.5]})
#ax1.bar(bin_centers, true_p, alpha = 0.5, color = 'red', width =bin_width,label = 'True distribution')
ax1.bar(bin_centers, reco_p, alpha = 0.5, color = 'blue', width = bin_width,label = 'Reco distribution')
ax1.errorbar(bin_centers, contents,
            yerr=errors,
            elinewidth=1,
            capsize=0,
            ls='None', marker='x', ms=5, color = 'black',
            label='Unfolded distribution')
# ax1.text(20, 450000, f'{ts} < {ts_st}')
# ax1.text(20, 400000, f'{numit} iter.')
ax1.set(xlabel='', ylabel='Counts')
ax1.legend()
#ax2.set_ylim(-5, 5)

#ax2.bar(bin_centers, pulls, alpha=0.5, color='grey', label='Pulls', width= bin_width)
#ax2.set_ylabel('Pulls')
#ax2.set_xlabel('$w$')
# ax2.set_xticks(bin_positions[::100]) 
#ax2.set_yticks([])
#ax2.axhline(0, color='black',linewidth=1) 
#ax2.axhline(-2.5, color='black',linewidth=1, ls = ':')
#ax2.axhline(2.5, color='black',linewidth=1, ls = ':')
plt.savefig(f'{figdir}/w_unfold_NOEFF_iter_{n_iter}.png', dpi=300)
plt.close()


#### plot normalized unfolded distribution 
plt.figure(figsize=(10, 8))
contents_normalized = contents / np.sum(contents)  # Normalize to unit area
errors_normalized = errors / np.sum(contents)  # Propagate errors for normalized distribution
bin_width = bin_centers[1] - bin_centers[0]
plt.errorbar(bin_centers, contents_normalized, xerr=bin_width/2, yerr=errors_normalized, fmt='none', color='blue', capsize=3)
plt.xlabel(r"$w$", fontsize=54)
plt.ylabel(r"$N^{\text{norm}}_{D^* \mu}$", fontsize=54)  # use LaTeX label
plt.tight_layout()
plt.savefig(f'{figdir}/w_unfold_NOEFF_normalized_iter_{n_iter}.png', dpi=300)
plt.close()

print(f"\n--- Unfolded Yields (Before Efficiency Correction) with nIter={n_iter} ---")
for w, y, e in zip(bin_centers, contents, errors):
#for w, y, e in zip(bin_centers, contents_normalized, errors_normalized):
    print(f"w = {w:.3f} | Yield = {y:.2f} +/- {e:.2f}")
#    print(f"w = {w:.3f} | Yield_norm = {y:.2f} +/- {e:.2f}")


print('Applying Efficiency Correction & Extracting Covariance...')

# 1. Load the Efficiency Histogram
eff_file_path = "/eos/lhcb/user/h/havva/q2Measurement-BtoDstMuNu/EffStudies/GenLevEff/Gen-lev-Eff-Combined.root"
eff_file = ROOT.TFile.Open(eff_file_path, "READ")
h_eff = eff_file.Get("gen_lev_eff_combined")

# Extract efficiency values (assuming the binning perfectly matches nbins=10, 1.0-1.5)
efficiencies = np.array([h_eff.GetBinContent(i) for i in range(1, nbins+1)])
eff_file.Close()



# 3. Extract the Covariance Matrix from RooUnfold
# Eunfold(1) returns the covariance matrix (ROOT.TMatrixD)
cov_matrix_root = unfold.Eunfold(ROOT.RooUnfold.kCovariance)
cov_matrix_unf = np.zeros((nbins, nbins))
cov_matrix_corr = np.zeros((nbins, nbins))

for i in range(nbins):
    for j in range(nbins):
        # TMatrixD is 0-indexed in Python
        cov_matrix_unf[i, j] = cov_matrix_root[i][j]
        # Scale the covariance by the efficiencies!
        cov_matrix_corr[i, j] = cov_matrix_unf[i, j] / (efficiencies[i] * efficiencies[j])


# Compute eff corrected unfolded histogram for plotting and saving to ROOT file
h_unfolded_corr = h_unfolded.Clone("h_unfolded_corr")
h_unfolded_corr.Divide(h_eff)  # Divide by efficiency bin-by-bin
retention_rate = 0.08  # integrated [Accepted / generated in simulation], adjust as needed
h_unfolded_corr.Scale(1/retention_rate)  # Scale by retention rate since we use different samples 1) from book-keeping 2)Gen-level production


#print the final yields after efficiency correction
print("\n--- Final Corrected and Unfolded Yields (Unblinded) ---")
unf_corr_values = np.array([h_unfolded_corr.GetBinContent(i) for i in range(1, nbins+1)])
unf_corr_errors = np.array([h_unfolded_corr.GetBinError(i) for i in range(1, nbins+1)])
normalised_unf_corr_values = unf_corr_values / np.sum(unf_corr_values)  # Normalize to unit area
normalised_unf_corr_errors = unf_corr_errors / np.sum(unf_corr_values)  # Propagate errors for normalized distribution
#for w, y, e, eff in zip(bin_centers, unf_corr_values, unf_corr_errors, efficiencies):
for w, y, e, eff in zip(bin_centers, normalised_unf_corr_values, normalised_unf_corr_errors, efficiencies):
    print(f"w = {w:.3f} | Eff = {eff*0.08*y:.5f} | Yield = {y:.5f} +/- {e:.5f}")

# Plot eff corrected unfolded distribution - blinded version (scale by blinding factor)
plt.figure(figsize=(10, 8))

blinding_array = [np.random.uniform(0.0, 1.0) for _ in range(nbins)]  # random numbers we cannot see

bin_centers = np.array([h_unfolded_corr.GetBinCenter(i) for i in range(1, nbins+1)])
contents_corr = np.array([h_unfolded_corr.GetBinContent(i) for i in range(1, nbins+1)])
contents_corr_blinded = contents_corr*blinding_array  # Apply blinding factor to the contents

errors_corr = np.array([h_unfolded_corr.GetBinError(i) for i in range(1, nbins+1)])
errors_corr_blinded = errors_corr*blinding_array  # Apply blinding factor to the errors

contents_normalized = contents_corr / np.sum(contents_corr)  # Normalize to unit area
errors_normalized = errors_corr / np.sum(contents_corr)  # Propagate errors for normalized distribution

contents_normalized_blinded = contents_corr_blinded / np.sum(contents_corr_blinded)  # Apply blinding factor to the normalized contents
errors_normalized_blinded = errors_corr_blinded / np.sum(contents_corr_blinded)  # Apply blinding factor to the normalized errors

w_width = (bin_centers[1] - bin_centers[0])*0.5

#plt.errorbar(bin_centers, contents_corr, xerr=w_width, yerr=errors_corr, fmt='none', color='blue', capsize=3)#, label='Efficiency-Corrected Unfolded')
plt.errorbar(bin_centers, contents_normalized_blinded, xerr=w_width, yerr=errors_normalized_blinded, fmt='none', color='blue', capsize=3)#, label='Efficiency-Corrected Unfolded')

plt.xlabel(r"$w$", fontsize=54)
#plt.ylabel(r"$N_{D^* \mu}$", fontsize=54)  # use LaTeX label
plt.ylabel(r"$N^{\text{norm}}_{D^* \mu}$", fontsize=54)  # use LaTeX label
#plt.title('Efficiency Corrected and Unfolded Distribution')
plt.tight_layout()
plt.savefig(f'{figdir}/w_unfold_eff_corrected_{nbins}_blinded.png', dpi=300)
plt.close()


# Print the final yields after efficiency correction
print("\n--- Final Corrected and Unfolded Yields (Blinded) ---")
for w, y,b, e, eff in zip(bin_centers, contents_corr_blinded, blinding_array, errors_corr_blinded, efficiencies):
    print(f"w = {w:.3f} | Eff = {eff*0.08*b:.3f} | Yield = {y:.2f} +/- {e:.2f}")


# ---------------------------------------------------------
# Save to a ROOT file 
# ---------------------------------------------------------
out_root_path = f"./unfolded_data_w_iteration_{n_iter}.root"
out_file = ROOT.TFile(out_root_path, "RECREATE")
h_unfolded.SetName("h_unfolded_data")
#h_unfolded.SetTitle("Unfolded w distribution")
h_unfolded.Write()
h_unfolded_corr.Write("h_unfolded_eff_corrected")
hresponse.Write("h_response_matrix")
# 3. Save the Efficiency-Corrected Covariance Matrix
# Create a 2D histogram with the same w binning as your 1D data
h_cov = ROOT.TH2D("h_cov_matrix", "Corrected Covariance Matrix;w;w", 
                  nbins, xmin, xmax, 
                  nbins, xmin, xmax)

# Fill the TH2D with values from your NumPy array
for i in range(nbins):
    for j in range(nbins):
        # Watch out for the indexing! 
        # NumPy is 0-indexed [i, j], but ROOT bins are 1-indexed (i+1, j+1)
        val = cov_matrix_corr[i, j]
        h_cov.SetBinContent(i + 1, j + 1, val)

h_cov.Write()
out_file.Close()
print(f"\nSaved unfolded TH1 to: {out_root_path}")