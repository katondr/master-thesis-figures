#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:51:08 2021

@author: tony
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as opt
import copy
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.sans-serif"] = "Times New Roman"
plt.rcParams['font.size'] = "9"

fig, axs = plt.subplots(nrows=5, ncols=3, sharex="row", sharey="row",
                           facecolor='white', figsize=(5.5, 8))

counter = 0
features=[
    GraphicFeature(start=0, end=20, strand=+1, color="tab:grey"),
    GraphicFeature(start=77, end=117, strand=-1, color="tab:red"),
    GraphicFeature(start=141, end=195, strand=+1, color="tab:blue"),
    GraphicFeature(start=204, end=215, strand=+1, color="tab:blue"),
    GraphicFeature(start=222, end=944, strand=+1, color="tab:olive"),
    GraphicFeature(start=953, end=1081, strand=+1, color="tab:red"),
    GraphicFeature(start=1223, end=1975, strand=-1, color="gold"),
    GraphicFeature(start=2021, end=2054, strand=-1, color="tab:blue"),
    GraphicFeature(start=2123, end=2177, strand=+1, color="tab:blue"),
    GraphicFeature(start=2194, end=2836, strand=+1, color="gold"),
    GraphicFeature(start=2858, end=2915, strand=+1, color="tab:red"),
    GraphicFeature(start=3011, end=3030, strand=-1, color="tab:grey"),
]
linear_features = copy.deepcopy(features)
#linear_features[0].label = "VF2"
#linear_features[0].fontdict["fontsize"] = 9
#linear_features[2].label = "$\it{luxI}$ promoter"
#linear_features[2].fontdict["fontsize"] = 9
#linear_features[3].label = "BBa_B0034"
#linear_features[3].fontdict["fontsize"] = 9
linear_features[4].label = "EYFP"
linear_features[4].fontdict["fontsize"] = 9
#linear_features[5].label = "terminator"
#linear_features[5].fontdict["fontsize"] = 9
linear_features[6].label = "$\it{luxR}$"
linear_features[6].fontdict["fontsize"] = 9
#linear_features[7].label = "$\it{luxR}$ promoter"
#linear_features[7].fontdict["fontsize"] = 9
#linear_features[8].label = "$\it{luxI}$ promoter"
#linear_features[8].fontdict["fontsize"] = 9
linear_features[9].label = "$\it{luxI}$"
linear_features[9].fontdict["fontsize"] = 9
#linear_features[11].label = "VR"
#linear_features[11].fontdict["fontsize"] = 9

gs = axs[0, 0].get_gridspec()
axs[0,0].remove()
axs[0,1].remove()
axs[0,2].remove()
record = GraphicRecord(sequence_length=3030, features=linear_features)
axbig = fig.add_subplot(gs[0,:])
record.plot(ax=axbig)
axbig.set_xticks(np.arange(0, 3200, 750))


plasmid_1_features = features.copy() + [
    GraphicFeature(start=3130, end=3718, strand=-1, color="tab:grey"),
    GraphicFeature(start=4016, end=4675, strand=-1, color="tab:green")

]
plasmid_2_features = features.copy() + [
    GraphicFeature(start=3477, end=4301, strand=+1, color="tab:grey"),
    GraphicFeature(start=4747, end=5409, strand=+1, color="tab:green")

]
plasmid_3_features = features.copy() + [
    GraphicFeature(start=3964, end=4914, strand=+1, color="tab:grey"),
    GraphicFeature(start=5230, end=5892, strand=+1, color="tab:green")

]
plasmid_1 = CircularGraphicRecord(sequence_length=4786, features=plasmid_1_features,
                                 top_position=3030/2)
plasmid_1.plot(ax=axs[1,2])

plasmid_2 = CircularGraphicRecord(sequence_length=5488, features=plasmid_2_features,
                                 top_position=3030/2)
plasmid_2.plot(ax=axs[1,1])

plasmid_3 = CircularGraphicRecord(sequence_length=5971, features=plasmid_3_features,
                                 top_position=3030/2)
plasmid_3.plot(ax=axs[1,0])

axs[1,0].text(-1.4, -2, "Cm$^R$", ha="center")
axs[1,0].text(0.4, -2.5, "pSC101 ori", ha="center")

axs[1,1].text(-1.4, -2, "Cm$^R$", ha="center")
axs[1,1].text(1, -2.4, "p15A ori", ha="center")

axs[1,2].text(-1, -2.4, "Cm$^R$", ha="center")
axs[1,2].text(1, -2.4, "pUC ori", ha="center")

data = pd.read_csv('../data/fl-210615.csv')
#time = data.t.unique().tolist()[:-1]
#time = ["0", "2", "4", "6", "8", "10", "12"]
time = [0, 2, 4, 6, 8, 10, 12]
label = data.x.unique().tolist()
abm = data["A"][data["t"] == "blc"].mean()
fbm = data["F"][data["t"] == "blc"].mean()

data["A"] -= abm
data["F"] -= fbm
data["F/A"] = data["F"]/data["A"]
data = data.loc[data["t"] != "blc"]
data["t"] = pd.to_numeric(data["t"])

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
    result.loc[index, "F/A_average"] = data["F/A"][
        (data["x"]==row["x"])&(data["t"]==row["t"])
    ].mean()
for index, row in result.iterrows():
    result.loc[index, "F/A_std"] = data["F/A"][
        (data["x"]==row["x"])&(data["t"]==row["t"])
    ].std()

x = np.arange(len(time))
width = 0.35

for column in np.arange(len(label)):
    axs[2, column].errorbar(x, 
                            result["A_average"][result["x"] == label[counter]],
                            yerr=result["A_std"][result["x"] == label[counter]],
                            ecolor="tab:grey",  color="tab:orange",
                            capsize=2, capthick=0.5, lw=1, mec="tab:grey",
                            marker=".", markersize=1)
    axs[2, column].spines["right"].set_visible(False)
    axs[2, column].spines["top"].set_visible(False)
    axs[2, column].set_xticks(x)
    axs[2, column].set_xticklabels(time)
    axs[2, column].set_ylim([0, 1])
    axs[2, column].set_yticks(np.arange(0, 1.2, 1))
    axs[2, column].yaxis.set_ticks_position('none')
    axs[2, column].xaxis.set_ticks_position('none')
    
    axs[3, column].errorbar(x,
                            result["F/A_average"][result["x"] == label[counter]],
                            yerr= result["F/A_std"][result["x"] == label[counter]],
                            ecolor="tab:grey",  color="tab:olive",
                            capsize=2, capthick=0.5, lw=1, mec="tab:gray",
                            marker=".", markersize=1)
    axs[3, column].spines["right"].set_visible(False)
    axs[3, column].spines["top"].set_visible(False)
    axs[3, column].set_xticks(x)
    axs[3, column].set_xticklabels(time)
    #axs[3, column].set_ylim([0, 120000000])
    #axs[3, column].set_yticks(np.arange(0, 120000000, 100000000))
    axs[3, column].yaxis.set_ticks_position('none')
    axs[3, column].xaxis.set_ticks_position('none')

    temp_data = result["A_average"][result["x"]==label[counter]]
    temp_data.reset_index(drop=True, inplace=True)
    masking = []
    for index in temp_data.index:
        #if index == temp_data.index[0]:
        #    if temp_data[index] < temp_data[index+1]:
        #        masking.append(index)      
        if index != temp_data.index[0]:
            if temp_data[index-1] < temp_data[index]:
                masking.append(index)
    bool_masking = bool()
    for item in masking:
        item_masking = (data["t"] == time[item]) & (data["x"] == label[counter])
        bool_masking = bool_masking | item_masking    
            
    x_axis = data["A"][(data["x"]==label[counter]) & bool_masking]
    y_axis = data["F/A"][(data["x"]==label[counter]) & bool_masking]
  
    axs[4, column].plot(x_axis, y_axis, "o", ms=2, color="tab:olive")
    axs[4, column].spines["right"].set_visible(False)
    axs[4, column].spines["top"].set_visible(False)
    #axs[4, column].set_ylim([0, 120000000])
    #axs[4, column].set_yticks(np.arange(0, 120000000, 100000000))
    axs[4, column].set_xlim([0, 1])
    axs[4, column].set_xticks(np.arange(0, 1.2, 1))
    axs[4, column].yaxis.set_ticks_position('none')
    axs[4, column].xaxis.set_ticks_position('none')
    def exp(x, a, b, y0):
        return (a * np.exp(b*x)) + y0
    #p0 = [70000, 9, min(y_axis)] # this is an mandatory initial guess
    popt, pcov = opt.curve_fit(exp, x_axis, y_axis, method="dogbox")
    x_fit = np.linspace(min(x_axis), max(x_axis), num= 300)
    y_fit = exp(x_fit, *popt)
    axs[4, column].plot(x_fit, y_fit, color="tab:blue", lw=1)
    print(popt)
    counter = counter + 1
    
axs[2,0].set_ylabel("A$_{600}$")
axs[3,0].set_ylabel("EYFP / A$_{600}$")
axs[4,0].set_ylabel("EYFP / A$_{600}$")
axs[2,0].set_xlabel("Hour")
axs[2,1].set_xlabel("Hour")
axs[2,2].set_xlabel("Hour")
axs[3,0].set_xlabel("Hour")
axs[3,1].set_xlabel("Hour")
axs[3,2].set_xlabel("Hour")
axs[4,0].set_xlabel("A$_{600}$")
axs[4,1].set_xlabel("A$_{600}$")
axs[4,2].set_xlabel("A$_{600}$")

axbig.text(-0.12, 0.5, "(a)", transform=axbig.transAxes, fontsize=9, va='top', ha='right')
axs[1,0].text(-0.8, 0.5, "(b)", transform=axs[1,0].transAxes, fontsize=9, va='top', ha='right')
axs[2,0].text(-0.4, 0.5, "(c)", transform=axs[2,0].transAxes, fontsize=9, va='top', ha='right')
axs[3,0].text(-0.4, 0.5, "(d)", transform=axs[3,0].transAxes, fontsize=9, va='top', ha='right')
axs[4,0].text(-0.4, 0.5, "(e)", transform=axs[4,0].transAxes, fontsize=9, va='top', ha='right')

fig.tight_layout()
#fig.savefig("../figure/thesis-34y-luxri.png", dpi=300, transparent=True)
#fig.savefig("../figure/thesis-34y-luxri.pdf")
fig.savefig("../figure/thesis-34y-luxri.png", dpi=300, bbox_inches='tight', transparent=True)