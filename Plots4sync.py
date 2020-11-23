# This code does the plotting for the paper "Synchronization between keyboard typing and neural oscillations"
# It also outputs dataframes used for statistical analyses

# Version 2 19/11/2020
# duprez.joan@gmail.com
## https://duprezjoan.wixsite.com/joanduprez ##

import numpy as np
import scipy.io
import scipy as sp
import matplotlib as mpl
from matplotlib import pyplot as plt
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
from plotly.offline import plot
import pandas as pd
import math
import seaborn as sns

outfold = 'yourfold'
## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
# RT results
# load the data

RTmat = scipy.io.loadmat(outfold + '/RTmat.mat')
RTdat = RTmat['RTmat']

arr = np.column_stack(
    list(map(np.ravel, np.meshgrid(*map(np.arange, RTdat.shape), indexing="ij"))) + [RTdat.ravel()])
RTmat2 = pd.DataFrame(arr, columns=['Subject', 'Condition', 'Precision', 'RT'])

RTmat2 = RTmat2.replace({'Precision': 0}, {'Precision': 'Correct trials'})
RTmat2 = RTmat2.replace({'Precision': 1}, {'Precision': 'Corrected errors'})
RTmat2 = RTmat2.replace({'Precision': 2}, {'Precision': 'Errors'})

RTmat2 = RTmat2.replace({'Condition': 0}, {'Condition': 'Words'})
RTmat2 = RTmat2.replace({'Condition': 1}, {'Condition': 'Pseudo-words'})
RTmat2 = RTmat2.replace({'Condition': 2}, {'Condition': 'Sentences'})
RTmat2 = RTmat2.replace({'Condition': 3}, {'Condition': 'Pseudo-sentences'})

# Export file for statistical analyses in R
RTmat2.to_csv(outfold + '/dataRT4stat.csv')

# Plot

# using seaborn and split violin plots

# create a copy of the data to put IKI in hz
RTmat3 = RTmat2.copy()
RTmat3 = RTmat3[RTmat3.Precision != 'Correct trials']
RTmat4 = RTmat3.copy()
ax.ylim = ([100, 1500])


colors = ['#1f77b4', '#9467bd', '#ff7f0e', '#d62728']

plt.subplot(1, 2, 1)

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

ax = [None] * (1 + 1)
fig = plt.figure(figsize=(12, 5))
gs = mpl.gridspec.GridSpec(nrows=1,
                           ncols=2,
                           figure=fig,
                           wspace=0.15, hspace=0.05
                           )

ax[0] = fig.add_subplot(gs[0, 0])
sns.violinplot(x="Condition", y="RT", data=RTmat2[RTmat2.Precision == 'Correct trials'][RTmat2.RT != 0])
sns.swarmplot(x='Condition', y='RT', data=RTmat2[RTmat2.Precision == 'Correct trials'][RTmat2.RT != 0], color="black",
              size=3, edgecolor='black', linewidth=0.5, alpha=0.5)

ax[0].collections[0].set_facecolor('tab:blue')
ax[0].collections[2].set_facecolor('tab:purple')
ax[0].collections[4].set_facecolor('tab:orange')
ax[0].collections[6].set_facecolor('tab:red')
ax[0].set_ylim([0, 2000])
ax[0].set_yticks([0, 500, 1000, 1500, 2000])
ax[0].set_ylabel('RT (ms)')
plt.xticks(fontsize=9)
ax[1] = fig.add_subplot(gs[0, 1])
sns.violinplot(x="Condition", y="RT", hue="Precision", data=RTmat3[RTmat3.RT != 0], split=True)
sns.swarmplot(x='Condition', y='RT', hue="Precision", data=RTmat3[RTmat3.RT != 0], dodge=True, color="white", size=3,
              edgecolor='black', linewidth=0.5, alpha=0.5)

ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[1].spines['left'].set_visible(False)
ax[1].axes.get_yaxis().set_visible(False)
ax[1].legend().set_visible(False)
plt.xticks(fontsize=9)
ax[1].collections[0].set_facecolor('tab:blue')
ax[1].collections[1].set_facecolor('lightsteelblue')
ax[1].collections[3].set_facecolor('tab:purple')
ax[1].collections[4].set_facecolor('plum')
ax[1].collections[6].set_facecolor('tab:orange')
ax[1].collections[7].set_facecolor('wheat')
ax[1].collections[9].set_facecolor('tab:red')
ax[1].collections[10].set_facecolor('salmon')

ax[0].set_xlabel('')
ax[1].set_xlabel('')
ax[0].set_title('Correct trials', fontweight='bold')
ax[1].set_title('Corrected errors (darker colors)\nand other errors (lighter colors)', fontweight='bold')

fig.savefig(outfold + '/RTviolin_split.png', format='png',
            dpi=1000)

## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
# IKI results

# for correct responses
# load the data
avgIKI_corr = scipy.io.loadmat(outfold + '/avgIKI_corr.mat')

test = pd.DataFrame(data=avgIKI_corr['avgIKI_corr'])

test2 = test
test2.columns = ["W", "pW", "S", "pS"]
test2.insert(2, "sub", np.linspace(1, 29, num=29))

test2 = pd.melt(test2, id_vars=['sub'],
                var_name='Condition', value_name='IKI')

test2 = test2.replace({'Condition': 'W'}, {'Condition': 'Words'})
test2 = test2.replace({'Condition': 'pW'}, {'Condition': 'Pseudo-words'})
test2 = test2.replace({'Condition': 'S'}, {'Condition': 'Sentences'})
test2 = test2.replace({'Condition': 'pS'}, {'Condition': 'Pseudo-sentences'})

# for backspace errors
avgIKI_BS = scipy.io.loadmat(outfold + '/avgIKI_BS.mat')

test = pd.DataFrame(data=avgIKI_BS['avgIKI_BS'])

test3 = test
test3.columns = ["W", "pW", "S", "pS"]
test3.insert(2, "sub", np.linspace(1, 29, num=29))

test3 = pd.melt(test3, id_vars=['sub'],
                var_name='Condition', value_name='IKI')
test3 = test3.replace({'Condition': 'W'}, {'Condition': 'Words'})
test3 = test3.replace({'Condition': 'pW'}, {'Condition': 'Pseudo-words'})
test3 = test3.replace({'Condition': 'S'}, {'Condition': 'Sentences'})
test3 = test3.replace({'Condition': 'pS'}, {'Condition': 'Pseudo-sentences'})

#  for other errors
avgIKI_err = scipy.io.loadmat(outfold + '/avgIKI_err.mat')

test = pd.DataFrame(data=avgIKI_err['avgIKI_err'])

test4 = test
test4.columns = ["W", "pW", "S", "pS"]
test4.insert(2, "sub", np.linspace(1, 29, num=29))

test4 = pd.melt(test4, id_vars=['sub'],
                var_name='Condition', value_name='IKI')
test4 = test4.replace({'Condition': 'W'}, {'Condition': 'Words'})
test4 = test4.replace({'Condition': 'pW'}, {'Condition': 'Pseudo-words'})
test4 = test4.replace({'Condition': 'S'}, {'Condition': 'Sentences'})
test4 = test4.replace({'Condition': 'pS'}, {'Condition': 'Pseudo-sentences'})

# format for stats on R

test2 = test2.assign(outcome='Correct trials')
test3 = test3.assign(outcome='Corrected errors')
test4 = test4.assign(outcome='Errors')

# Concatenate pd frames
frames = [test2, test3, test4]
dataIKI = pd.concat(frames)
dataIKI.to_csv(outfold + '/dataIKI4stat.csv')

# plot using seaborn and split violin plots

# create a copy of the data to put IKI in hz
dataIKI2 = dataIKI.copy()
dataIKI2 = dataIKI2[dataIKI2.outcome != 'Correct trials']
dataIKI3 = dataIKI2.copy()
ax.ylim = ([100, 300])
for triali in np.linspace(0, 215, 216, dtype='int64'):
    dataIKI3.IKI.iloc[triali] = 1000 / dataIKI3.IKI.iloc[triali]

colors = ['#1f77b4', '#9467bd', '#ff7f0e', '#d62728']

plt.subplot(1, 2, 1)

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

ax = [None] * (1 + 1)
fig = plt.figure(figsize=(12, 5))
gs = mpl.gridspec.GridSpec(nrows=1,
                           ncols=2,
                           figure=fig,
                           wspace=0.15, hspace=0.05
                           )

ax[0] = fig.add_subplot(gs[0, 0])
sns.violinplot(x="Condition", y="IKI", data=dataIKI[dataIKI.outcome == 'Correct trials'])
sns.swarmplot(x='Condition', y='IKI', data=dataIKI[dataIKI.outcome == 'Correct trials'], color="black", size=3,
              edgecolor='black', linewidth=0.5, alpha=0.5)

ax[0].collections[0].set_facecolor('tab:blue')
ax[0].collections[2].set_facecolor('tab:purple')
ax[0].collections[4].set_facecolor('tab:orange')
ax[0].collections[6].set_facecolor('tab:red')
ax[0].set_ylim([50, 350])
ax[0].set_yticks([100, 200, 300])
ax[0].set_ylabel('IKI (ms)')
plt.xticks(fontsize=9)

ax[1] = fig.add_subplot(gs[0, 1])
sns.violinplot(x="Condition", y="IKI", hue="outcome", data=dataIKI2, split=True)
sns.swarmplot(x='Condition', y='IKI', hue="outcome", data=dataIKI2, dodge=True, color="white", size=3,
              edgecolor='black', linewidth=0.5, alpha=0.5)
plt.xticks(fontsize=9)

ax2 = plt.twinx()
ax2.plot()
ax2.set_ylim(ax[0].get_ylim())
ax2.set_ylim([50, 350])
ax2.set_yticks([100, 200, 300])
ax2.set_yticklabels(['10', '5', '3.3'])
ax2.set_ylabel('Frequency (Hz)')

ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[1].spines['left'].set_visible(False)
ax[1].axes.get_yaxis().set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax[1].legend().set_visible(False)

ax[1].collections[0].set_facecolor('tab:blue')
ax[1].collections[1].set_facecolor('lightsteelblue')
ax[1].collections[3].set_facecolor('tab:purple')
ax[1].collections[4].set_facecolor('plum')
ax[1].collections[6].set_facecolor('tab:orange')
ax[1].collections[7].set_facecolor('wheat')
ax[1].collections[9].set_facecolor('tab:red')
ax[1].collections[10].set_facecolor('salmon')

ax[0].set_xlabel('')
ax[1].set_xlabel('')
ax[0].set_title('Correct trials', fontweight='bold')
ax[1].set_title('Corrected errors (darker colors)\nand other errors (lighter colors)', fontweight='bold')

fig.savefig(outfold + '/IKIviolin_split.png', format='png',
            dpi=1000)

## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
# KDE analysis results

# load KDE data
kderes = scipy.io.loadmat(outfold + '/kde_res4py.mat')
hz = np.linspace(0, 100, 1000)

condcorr = ['hz_correctW', 'hz_correctpseudoW', 'hz_correctS', 'hz_correctpseudoS']
kdecorr = np.ndarray(shape=(27, 1000, 4))  # use kderes['hz_BSS'].shape[0] to soft code the size of ndarray dimensions
for condi in np.arange(0, 4, 1):
    kdecorr[:, :, condi] = kderes[condcorr[condi]]

# Plot grand average result

condall = ['hz_correctW', 'hz_correctpseudoW', 'hz_correctS', 'hz_correctpseudoS', 'hz_BSW', 'hz_BSpseudoW', 'hz_BSS',
           'hz_BSpseudoS',
           'hz_errW', 'hz_errpseudoW', 'hz_errS', 'hz_errpseudoS']
kdeall = np.ndarray(shape=(27, 1000, 12))  # use kderes['hz_BSS'].shape[0] to soft code the size of ndarray dimensions
for condi in np.arange(0, 12, 1):
    kdeall[:, :, condi] = kderes[condall[condi]]

avgkdeall = np.mean(kdeall, 2)

# Get averaged peak frequency for all subs

peakfreqKDE = np.empty(27)
for subi in np.linspace(0, 26, 27, dtype='int64'):
            peakfreqKDE[subi,] = hz[avgkdeall[subi,].argmax()]

kdecorravg = np.mean(kdecorr, axis=2)

condbs = ['hz_BSW', 'hz_BSpseudoW', 'hz_BSS', 'hz_BSpseudoS']
kdebs = np.ndarray(shape=(27, 1000, 4))  # use kderes['hz_BSS'].shape[0] to soft code the size of ndarray dimensions
for condi in np.arange(0, 4, 1):
    kdebs[:, :, condi] = kderes[condbs[condi]]

kdebsavg = np.mean(kdebs, axis=2)

conderr = ['hz_errW', 'hz_errpseudoW', 'hz_errS', 'hz_errpseudoS']
kdeerr = np.ndarray(shape=(27, 1000, 4))  # use kderes['hz_BSS'].shape[0] to soft code the size of ndarray dimensions
for condi in np.arange(0, 4, 1):
    kdeerr[:, :, condi] = kderes[conderr[condi]]

kdeerravg = np.mean(kdeerr, axis=2)

# Peak frequency KDE analysis

peakf_c = scipy.io.loadmat(outfold + '/peakFreq_corr.mat')
pfreq = pd.DataFrame(data=peakf_c['pfcorr'])
pfreq.columns = ["W", "pW", "S", "pS"]
pfreq.insert(2, "sub", np.linspace(1, len(pfreq), num=len(pfreq)))
# Tidy the dataframe for plotting using melt of pandas
pfreq = pd.melt(pfreq, id_vars=['sub'],
                var_name='Condition', value_name='peak frequency')


# for corrected errors
peakf_bs = scipy.io.loadmat(outfold + '/peakFreq_BS.mat')
pfreq_bs = pd.DataFrame(data=peakf_bs['pfBS'])
pfreq_bs.columns = ["W", "pW", "S", "pS"]
pfreq_bs.insert(2, "sub", np.linspace(1, len(pfreq_bs), num=len(pfreq_bs)))
# Tidy the dataframe for plotting using melt of pandas
pfreq_bs = pd.melt(pfreq_bs, id_vars=['sub'],
                   var_name='Condition', value_name='peak frequency')


# for other errors
peakf_err = scipy.io.loadmat(outfold + '/peakFreq_err.mat')
pfreq_err = pd.DataFrame(data=peakf_err['pferr'])
pfreq_err.columns = ["W", "pW", "S", "pS"]
pfreq_err.insert(2, "sub", np.linspace(1, len(pfreq_err), num=len(pfreq_err)))
# Tidy the dataframe for plotting using melt of pandas
pfreq_err = pd.melt(pfreq_err, id_vars=['sub'],
                    var_name='Condition', value_name='peak frequency')

# format 4 stats on R

pfreq = pfreq.assign(outcome='Correct trials')
pfreq_bs = pfreq_bs.assign(outcome='Corrected errors')
pfreq_err = pfreq_err.assign(outcome='Errors')

# Concatenate pd frames
frames = [pfreq, pfreq_bs, pfreq_err]
dataKDE = pd.concat(frames)
dataKDE = dataKDE.replace({'Condition': 'W'}, {'Condition': 'Words'})
dataKDE = dataKDE.replace({'Condition': 'pW'}, {'Condition': 'Pseudo-words'})
dataKDE = dataKDE.replace({'Condition': 'S'}, {'Condition': 'Sentences'})
dataKDE = dataKDE.replace({'Condition': 'pS'}, {'Condition': 'Pseudo-sentences'})
dataKDE = dataKDE.rename(columns={'peak frequency': 'freq'})
dataKDE.to_csv(outfold + 'Dropbox/SINS/Keyboard/Keysync/dataKDE4stat.csv')

# Plot
ax = [None] * (1 + 1)
fig = plt.figure(figsize=(12, 5))

gs = mpl.gridspec.GridSpec(nrows=3,
                           ncols=2,
                           figure=fig,
                           wspace=0.15, hspace=1.15
                           )
fig.add_subplot(gs[0, 0])
plt.imshow(avgkdeall[:, :250], aspect='auto', cmap='inferno', extent=[0, 25, 0, 27])
plt.title('KDE spectrum')
plt.xlabel('Hz')
plt.ylabel('Subjects')
plt.text(-3, 30, 'A', fontsize=16, fontweight='bold')

ax[0] = fig.add_subplot(gs[1:, 0])
sns.violinplot(x="Condition", y="freq", data=dataKDE[dataKDE.outcome == 'Correct trials'])
sns.swarmplot(x='Condition', y='freq', data=dataKDE[dataKDE.outcome == 'Correct trials'], color="black", size=3,
              edgecolor='black', linewidth=0.5, alpha=0.5)

ax[0].collections[0].set_facecolor('tab:blue')
ax[0].collections[2].set_facecolor('tab:purple')
ax[0].collections[4].set_facecolor('tab:orange')
ax[0].collections[6].set_facecolor('tab:red')
ax[0].set_ylim([1, 15])
ax[0].set_yticks([1, 5, 10, 15])
ax[0].set_ylabel('Frequency (Hz)')
plt.xticks(fontsize=9)

ax[1] = fig.add_subplot(gs[1:, 1])

sns.violinplot(x="Condition", y="freq", hue="outcome", data=dataKDE[dataKDE.outcome != 'Correct trials'], split=True)
sns.swarmplot(x='Condition', y='freq', hue="outcome", data=dataKDE[dataKDE.outcome != 'Correct trials'], dodge=True,
              color="white", size=3, edgecolor='black', linewidth=0.5, alpha=0.5)
ax[1].set_ylim([1, 15])
ax[1].set_yticks([1, 5, 10, 15])
plt.xticks(fontsize=9)

ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[1].spines['left'].set_visible(False)
ax[1].axes.get_yaxis().set_visible(False)
ax[1].legend().set_visible(False)
ax[1].collections[0].set_facecolor('tab:blue')
ax[1].collections[1].set_facecolor('lightsteelblue')
ax[1].collections[3].set_facecolor('tab:purple')
ax[1].collections[4].set_facecolor('plum')
ax[1].collections[6].set_facecolor('tab:orange')
ax[1].collections[7].set_facecolor('wheat')
ax[1].collections[9].set_facecolor('tab:red')
ax[1].collections[10].set_facecolor('salmon')

ax[0].set_xlabel('')
ax[1].set_xlabel('')
ax[0].set_title('Correct trials', fontweight='bold')
ax[1].set_title('Corrected errors (darker colors)\nand other errors (lighter colors)', fontweight='bold')
ax[0].text(-0.95, 16, 'B', fontsize=16, fontweight='bold')

fig.savefig(outfold + '/KDEspect&violin_split.png',
            format='png', dpi=1000)

## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
# Topomaps of the GED components
# Plot the subject average topomaps of each frequency's component

# Get the GED spatial filter data
avgged = scipy.io.loadmat(outfold + '/avg_mapsTH.mat')
avg_mapsTH = avgged['avg_mapsTH']  # sub*freuency*channels

import mne

# Get the chanlocs by loading an .set file. Here is a special one from which electrodes 61 to 64 were removed previously
# in the EEG.chanlocs in matlab because MNE plots chan 61 around Cz. Don't know why. File can be found here : https://github.com/jduprez
raw = mne.io.read_raw_eeglab(
    outfold + 'test_chanloc2.set')  # only use to get chanlocs

# Define frequencies used for GED
frex = np.linspace(1, 15, 29)

# loop through GED frequencies to get the topomap and save as pdf
ploti = 1
for freqi in np.linspace(6, 28, 12, dtype='int64'):
    plt.subplot(3, 4, ploti)
    avg = np.nanmean(avg_mapsTH[:, freqi, :], axis=0)
    mne.viz.plot_topomap(avg.take(np.linspace(0, 59, 60, dtype='int64')), raw.info, sensors=False,
                         show=True, outlines='skirt', contours=0)
    plt.title(str(int(frex[freqi])) + ' Hz')
    ploti = ploti + 1

plt.subplots_adjust(hspace=0.5)
plt.savefig(outfold + 'Dropbox/SINS/Keyboard/GEDtopomaps_withoutsensors.png', format='png', dpi=1000)

# Plot with eigen values

eigval = scipy.io.loadmat(outfold + 'Dropbox/SINS/Keyboard/Keysync/mat4eval.mat')
mat4eval = eigval['mat4eval']

plt.plot(frex[4:], np.mean(mat4eval[:, 4:], 0))
plt.xlim(4, 15)
plt.xlabel('Frequency (Hz)', Fontsize=16)
plt.ylabel('Average highest eigenvalue', Fontsize=16)
plt.yticks([0.05, 0.10, 0.15, 0.20])

# Plot
fig = plt.figure(figsize=(6, 8))
ax = [None] * (1 + 1)
gs = mpl.gridspec.GridSpec(nrows=4,
                           ncols=5,
                           figure=fig,
                           left=0.05,
                           right=0.75,
                           wspace=0.1, hspace=0.75
                           )
ax[0] = fig.add_subplot(gs[0, 1:4])

plt.plot(frex[4:], np.mean(mat4eval[:, 4:], 0))
plt.xlim(4, 15)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Average \nhighest eigenvalue')
plt.yticks([0.10, 0.20])

freqi = 6
for row in range(3):
    for col in range(4):
        ax = fig.add_subplot(gs[row + 1, col])
        avg = np.nanmean(avg_mapsTH[:, freqi, :], axis=0)
        # fig = plt.figure()
        mne.viz.plot_topomap(avg.take(np.linspace(0, 59, 60, dtype='int64')), raw.info, sensors=False,
                             show=True, outlines='skirt', contours=0)
        plt.title(str(int(frex[freqi])) + ' Hz')

        freqi = freqi + 2
plt.text(-5.5, 7.25, 'A', fontsize=16, fontweight='bold')
plt.text(-5.5, 5.5, 'B', fontsize=16, fontweight='bold')

#ax = fig.add_subplot(gs[row + 1, 4])
fig, ax = plt.subplots(figsize=(1, 3))
fig.subplots_adjust(right=0.25)

cmap = mpl.cm.RdBu_r
norm = mpl.colors.Normalize(vmin=-3, vmax=3)

cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='vertical', ticks=[-3, 3])
cb1.set_label('a.u')
cb1.ax.set_yticklabels(['-', '+'])
fig.show()
plt.savefig(outfold + '/Eval&GEDtopomaps_withoutsensors.png', format='png', dpi=1000)


## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
# Plot condition specific phase clustering across the different GED frequencies

# Get the clustering data
# The data is a 3 (precision: correct, BSerr, err) * 25 (frequencies from 3 to 15Hz) * 4 (condition: W, pW...) * 30 (sub) * 2 (measure: itcp, itpcZ) matrix
clust = scipy.io.loadmat(outfold + '/clust_mat.mat')
clust_mat = clust['clust_mat']

# format 4 stats on R
# create empty dataframe
# get only ITPC

arr = np.column_stack(
    list(map(np.ravel, np.meshgrid(*map(np.arange, clust_mat.shape), indexing="ij"))) + [clust_mat.ravel()])
clust4stat = pd.DataFrame(arr, columns=['Precision', 'Frequency', 'Condition', 'Subject', 'ITPCz'])

clust4stat = clust4stat.replace({'Precision': 0}, {'Precision': 'Correct trials'})
clust4stat = clust4stat.replace({'Precision': 1}, {'Precision': 'Corrected errors'})
clust4stat = clust4stat.replace({'Precision': 2}, {'Precision': 'Errors'})

clust4stat = clust4stat.replace({'Condition': 0}, {'Condition': 'Words'})
clust4stat = clust4stat.replace({'Condition': 1}, {'Condition': 'Pseudowords'})
clust4stat = clust4stat.replace({'Condition': 2}, {'Condition': 'Sentences'})
clust4stat = clust4stat.replace({'Condition': 3}, {'Condition': 'Pseudosentences'})

clust4stat.to_csv(outfold + '/dataCLUST4stat_perm.csv')


# Plot

ax = [None] * (1 + 1)
fig = plt.figure(figsize=(10, 6))
gs = mpl.gridspec.GridSpec(nrows=3,
                           ncols=10,
                           figure=fig,
                           wspace=0.75, hspace=0.5
                           )

ax[0] = fig.add_subplot(gs[0, 0:4])

zclust_avg = np.transpose(np.mean(np.mean(np.mean(clust_mat[:, 6:29, :, :], 0), 1), 1))
#zclust_sem = np.transpose(scipy.stats.sem(np.mean(np.mean(clust_mat[:, 6:29, :, :], 0), 1), 1))
x = np.linspace(4, 15, 23)
plt.plot(x, np.transpose(np.mean(np.mean(clust_mat[0, 6:29, :, :], 2), 1)), label='Correct responses', linewidth=1)
plt.plot(x, np.transpose(np.mean(np.mean(clust_mat[1, 6:29, :, :], 2), 1)), label='Corrected errors', linewidth=1)
plt.plot(x, np.transpose(np.mean(np.mean(clust_mat[2, 6:29, :, :], 2), 1)), label='Other errors', linewidth=1)
plt.ylim(-0.3, 0.3)
plt.xticks([4, 6, 8, 10, 12, 14])
plt.xlabel('Frequency (Hz)')
plt.plot(x, np.transpose(np.mean(np.mean(np.mean(clust_mat[:, 6:29, :, :], 0), 1), 1)), 'k', linewidth=2, label='Average')
plt.xticks([4, 6, 8, 10, 12, 14])
plt.xlim([4, 15])
plt.ylabel('Z-scored synchronization')
leg = plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2,
            borderaxespad=0, frameon=False, fontsize=8)
ax[0].xaxis.set_ticklabels([])

ax[1] = fig.add_subplot(gs[1:3, 0:4])
tempclust = np.transpose(np.mean(np.mean(clust_mat[:, 6:29, :, :], 0), 1))

# Plot subject-specific spectrum of clustering
temp = tempclust[~np.all(tempclust == 0, axis=1)]

# Get only pixels of significant sync at the subject level

sigfreq = scipy.io.loadmat(bpath + '/sigfreq.mat')
sig4plot = sigfreq['sigfreq']
sig4plot = np.delete(sig4plot, np.s_[0:6], 1)
sig4plot = np.delete(sig4plot, np.s_[0, 3, 18], 0)

temp[sig4plot > 0.05] = 'nan'

peakfreq = np.empty([27])
frex = np.linspace(4, 15, 23)

        for subi in np.linspace(0, 26, 27, dtype='int64'):
            peakfreq[subi,] = frex[np.nanargmax(temp[subi,])]

# Sort subject according to peak freq
sorted_clust = temp[np.argsort(peakfreq),:]

mpl.pyplot.imshow(sorted_clust, cmap = 'jet', extent=[4,15,27,0], aspect='auto')
plt.clim(-1, 1)
plt.xlim([4, 15])
plt.ylabel('Subjects')
plt.xlabel('Frequency (Hz)')
cbaxes = fig.add_axes([0.425, 0.18, 0.01, 0.35])
plt.colorbar(cax = cbaxes)
plt.clim(-1, 1)
cbaxes.set_ylabel('Z-scored synchronization')

# Scatter plot with peak KDE frequency, need peakKDE !
# Recompute peakfreq because we don't need to sort only significant clustering values for correlation with behavior
temp = tempclust[~np.all(tempclust == 0, axis=1)]

peakfreq = np.empty([27])
frex = np.linspace(4, 15, 23)

        for subi in np.linspace(0, 26, 27, dtype='int64'):
            peakfreq[subi,] = frex[np.nanargmax(temp[subi,])]


ax[1] = fig.add_subplot(gs[1:3, 6:10])
g1 = plt.scatter(peakfreq,peakfreqKDE, marker= 'o', facecolors='none', edgecolors= 'C0', label="KDE")
plt.xlabel('Peak synchronization frequency (Hz)')
plt.ylabel('Behavioral peak frequency (Hz)')

plt.plot(np.unique(peakfreq), np.poly1d(np.polyfit(peakfreq, peakfreqKDE, 1))(np.unique(peakfreq)))

## Same with IKI defined frequency, needs dataIKI !
temp=dataIKI.groupby('sub').mean()
peakfreqIKI = temp.values
peakfreqIKI = 1000/peakfreqIKI

g2 = plt.scatter(peakfreq,peakfreqIKI, marker= 'o', facecolors='none', edgecolors= 'red', label='IKI')
plt.plot(np.unique(peakfreq), np.poly1d(np.squeeze(np.polyfit(peakfreq, peakfreqIKI, 1)))(np.unique(peakfreq)), color='red')
plt.legend([g1, g2], ['KDE, \u03C1 = 0.49, p = 0.008', 'IKI, \u03C1 = 0.35, p = 0.069'], frameon=False, bbox_to_anchor= (-0.1,0.18,1,1))
plt.xlim(3.5,15.5)

plt.text(-20, 10.1, 'A', fontsize=16, fontweight='bold')
plt.text(-20, 7.5, 'B', fontsize=16, fontweight='bold')
plt.text(-0.25, 7.5, 'C', fontsize=16, fontweight='bold')

plt.savefig(outfold + '/clustering-sub&gplevel.png', format='png', dpi=1000)

# END