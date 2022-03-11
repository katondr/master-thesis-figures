#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:47:15 2021

@author: tony
"""
#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy
from matplotlib.ticker import ScalarFormatter
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.sans-serif"] = "Times New Roman"
plt.rcParams['font.size'] = "9"
#%%
features=[
    GraphicFeature(start=0, end=20, strand=+1, color="tab:grey"),
    GraphicFeature(start=77, end=117, strand=-1, color="tab:red"),
    GraphicFeature(start=241, end=307, strand=+1, color="tab:blue"),
    GraphicFeature(start=349, end=360, strand=+1, color="tab:blue"),
    GraphicFeature(start=367, end=1044, strand=+1, color="tab:red"),
    GraphicFeature(start=1089, end=1160, strand=+1, color="tab:red"),
    GraphicFeature(start=1231, end=1288, strand=+1, color="tab:red"),
    GraphicFeature(start=1364, end=1383, strand=-1, color="tab:grey")
]
linear_features = copy.deepcopy(features)
linear_features[0].label = "VF2"
linear_features[0].fontdict["fontsize"] = 9
linear_features[3].label = "B0034"
linear_features[3].fontdict["fontsize"] = 9
linear_features[2].label = "lac promoter"
linear_features[2].fontdict["fontsize"] = 9
linear_features[4].label = "mRFP"
linear_features[4].fontdict["fontsize"] = 9
linear_features[5].label = "nrB T1"
linear_features[5].fontdict["fontsize"] = 9
linear_features[7].label = "VR"
linear_features[7].fontdict["fontsize"] = 9
record = GraphicRecord(sequence_length=1383, features=linear_features)
plasmid_1_features = features.copy() + [
    GraphicFeature(start=1483, end=2071, strand=-1, color="tab:grey"),
    GraphicFeature(start=2369, end=3028, strand=-1, color="tab:green")

]
plasmid_2_features = features.copy() + [
    GraphicFeature(start=1813, end=2357, strand=+1, color="tab:grey"),
    GraphicFeature(start=3066, end=3725, strand=+1, color="tab:green")

]
plasmid_3_features = features.copy() + [
    GraphicFeature(start=2283, end=3233, strand=+1, color="tab:grey"),
    GraphicFeature(start=3549, end=4211, strand=+1, color="tab:green")

]
plasmid_1 = CircularGraphicRecord(sequence_length=3139, features=plasmid_1_features,
                                 top_position=1383/2)
plasmid_2 = CircularGraphicRecord(sequence_length=3807, features=plasmid_2_features,
                                 top_position=1383/2)
plasmid_3 = CircularGraphicRecord(sequence_length=4290, features=plasmid_3_features,
                                 top_position=1383/2)
#%%
data = pd.read_csv('../data/fl-210411.csv')

time = data.t.unique().tolist()[1:]
label = data.x.unique().tolist()

abm = data["A"][data["t"] == "blc"].mean()
fbm_r = data["F"][
    (data["t"] == "blc") & 
    ((data["x"] == "1C") | (data["x"] == "3C") | (data["x"] == "4C"))
].mean()
fbm_y = data["F"][
    (data["t"] == "blc") & 
    ((data["x"] != "1C") | (data["x"] != "3C") | (data["x"] != "4C"))
].mean()

data["A"] -= abm
data.loc[(data["x"] == "1C") | (data["x"] == "3C") | (data["x"] == "4C"), ("F")]  -= fbm_r
data.loc[(data["x"] != "1C") | (data["x"] != "3C") | (data["x"] != "4C"), ("F")] -= fbm_y
data["F/A"] = data["F"]/data["A"]

label_r = label[:3]
label_r

result = pd.DataFrame({"x":[], "t":[], "A_average":[], "A_std":[],
                       "F/A_average":[], "F/A_std":[]})

temp = pd.DataFrame({"x":[], "t":[]})
for item in label:
    temp["t"] = time
    temp["x"] = [item]*7
    result = pd.concat([result, temp], ignore_index=True)

for index, row in result.iterrows():
    result.loc[index, "A_average"] = data["A"][
        (data["x"]==row["x"])&(data["t"]==row["t"])
    ].mean()
for index, row in result.iterrows():
    result.loc[index, "A_std"] = data["A"][
        (data["x"]==row["x"])&(data["t"]==row["t"])
    ].std()
for index, row in result.iterrows():
    result.loc[index, "F/A_average"] = data["F"][
        (data["x"]==row["x"])&(data["t"]==row["t"])
    ].mean()
for index, row in result.iterrows():
    result.loc[index, "F/A_std"] = data["F"][
        (data["x"]==row["x"])&(data["t"]==row["t"])
    ].std()

#%%
fig, axs = plt.subplots(nrows=5, ncols=3, sharex=False, sharey="row",
                           facecolor='white', figsize=(6, 9))
gs = axs[0, 0].get_gridspec()
axs[0,0].remove()
axs[0,1].remove()
axs[0,2].remove()
plasmid_1.plot(ax=axs[1,2])
plasmid_2.plot(ax=axs[1,1])
plasmid_3.plot(ax=axs[1,0])
axbig = fig.add_subplot(gs[0,:])
record.plot(ax=axbig)
counter = 0
x = np.arange(len(time))
width = 0.35
for column in np.arange(len(label_r)):
    axs[1, column].errorbar(x, 
                            result["A_average"][result["x"] == label_r[counter]],
                            yerr=result["A_std"][result["x"] == label_r[counter]],
                            ecolor="tab:grey",  color="tab:orange",
                            capsize=2, capthick=0.5, lw=1, mec="tab:grey",
                            marker=".", markersize=1)
    axs[1, column].spines["right"].set_visible(False)
    axs[1, column].spines["top"].set_visible(False)
    axs[1, column].set_xticks(x)
    axs[1, column].set_xticklabels(time)
    axs[1, column].set_ylim([0, 1])
    axs[1, column].set_yticks(np.arange(0, 1.2, 1))
    axs[1, column].yaxis.set_ticks_position('none')
    axs[1, column].xaxis.set_ticks_position('none')

    axs[2, column].errorbar(x,
                            result["F/A_average"][result["x"] == label_r[counter]],
                            yerr= result["F/A_std"][result["x"] == label_r[counter]],
                            ecolor="tab:grey",  color="tab:red",
                            capsize=2, capthick=0.5, lw=1, mec="tab:gray",
                            marker=".", markersize=1) 
    axs[2, column].spines["right"].set_visible(False)
    axs[2, column].spines["top"].set_visible(False)
    axs[2, column].set_xticks(x)
    axs[2, column].set_xticklabels(time)
    axs[2, column].set_ylim([0, 10000000])
    axs[2, column].set_yticks(np.arange(0, 12000000, 10000000))
    axs[2, column].yaxis.set_ticks_position('none')
    axs[2, column].xaxis.set_ticks_position('none')

    a, b = np.polyfit(data["A"][data["x"]==label_r[counter]], data["F/A"][data["x"]==label_r[counter]],deg=1)
    y_est = a * data["A"][data["x"]==label_r[counter]] + b
    axs[4, column].plot(data["A"][data["x"]==label_r[counter]],
                        data["F/A"][data["x"]==label_r[counter]], 
                        "o", ms=2, color="tab:red")
    axs[3, column].plot(data["A"][data["x"]==label_r[counter]], y_est,
                        color="tab:blue", lw=1)
    axs[3, column].spines["right"].set_visible(False)
    axs[3, column].spines["top"].set_visible(False)
    axs[3, column].set_ylim([0, 10000000])
    axs[3, column].set_yticks(np.arange(0, 12000000, 10000000))
    axs[3, column].set_xlim([0, 1])
    axs[3, column].set_xticks(np.arange(0, 1.2, 1))
    axs[3, column].yaxis.set_ticks_position('none')
    axs[3, column].xaxis.set_ticks_position('none')
    counter = counter + 1

#axs[1,0].text(-1.5, -1.4, "Cm$^R$", ha="center")
#axs[1,0].text(0.2, -2.5, "pSC101 ori", ha="center")
#axs[1,1].text(-1.5, -1.6, "Cm$^R$", ha="center")
#axs[1,1].text(1.5, -2, "p15A ori", ha="center")
#axs[1,2].text(-1.2, -2.1, "Cm$^R$", ha="center")
#axs[1,2].text(1.5, -2, "pUC ori", ha="center")
axs[1,0].set_ylabel("A$_{600}$")
axs[2,0].set_ylabel("mRFP / A$_{600}$")
axs[3,0].set_ylabel("mRFP / A$_{600}$")
axs[1,0].set_xlabel("Hour")
axs[1,1].set_xlabel("Hour")
axs[1,2].set_xlabel("Hour")
axs[2,0].set_xlabel("Hour")
axs[2,1].set_xlabel("Hour")
axs[2,2].set_xlabel("Hour")
axs[3,0].set_xlabel("A$_{600}$")
axs[4,1].set_xlabel("A$_{600}$")
axs[4,2].set_xlabel("A$_{600}$")
axbig.text(-0.11, 0.5, "(a)", transform=axbig.transAxes, fontsize=9, va='top', ha='right')
axs[1,0].text(-0.7, 0.5, "(b)", transform=axs[1,0].transAxes, fontsize=9, va='top', ha='right')
axs[2,0].text(-0.4, 0.5, "(c)", transform=axs[2,0].transAxes, fontsize=9, va='top', ha='right')
axs[3,0].text(-0.4, 0.5, "(d)", transform=axs[3,0].transAxes, fontsize=9, va='top', ha='right')
#axs[4,0].text(-0.4, 0.5, "(e)", transform=axs[4,0].transAxes, fontsize=9, va='top', ha='right')

fig.tight_layout()
#fig.savefig("../figure/thesis-mrfp.pdf")