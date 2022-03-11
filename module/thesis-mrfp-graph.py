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
fig, axs = plt.subplots(nrows=3, ncols=3, sharex=False, sharey="row",
                           facecolor='white', figsize=(5, 5.5))
counter = 0
x = np.arange(len(time))
width = 0.35
row = 0
for column in np.arange(len(label_r)):
    axs[0, column].errorbar(x, 
                            result["A_average"][result["x"] == label_r[counter]],
                            yerr=result["A_std"][result["x"] == label_r[counter]],
                            ecolor="tab:grey",  color="tab:orange",
                            capsize=2, capthick=0.5, lw=1, mec="tab:grey",
                            marker=".", markersize=1)
    axs[0, column].spines["right"].set_visible(False)
    axs[0, column].spines["top"].set_visible(False)
    axs[0, column].set_xticks(x)
    axs[0, column].set_xticklabels(time)
    axs[0, column].set_ylim([0, 1])
    axs[0, column].set_yticks(np.arange(0, 1.2, 1))
    #axs[0, column].yaxis.set_ticks_position('none')
    #axs[0, column].xaxis.set_ticks_position('none')

    axs[1, column].errorbar(x,
                            result["F/A_average"][result["x"] == label_r[counter]],
                            yerr= result["F/A_std"][result["x"] == label_r[counter]],
                            ecolor="tab:grey",  color="tab:red",
                            capsize=2, capthick=0.5, lw=1, mec="tab:gray",
                            marker=".", markersize=1) 
    axs[1, column].spines["right"].set_visible(False)
    axs[1, column].spines["top"].set_visible(False)
    axs[1, column].set_xticks(x)
    axs[1, column].set_xticklabels(time)
    axs[1, column].set_ylim([0, 10000000])
    axs[1, column].set_yticks(np.arange(0, 12000000, 10000000))
    #axs[1, column].yaxis.set_ticks_position('none')
    #axs[1, column].xaxis.set_ticks_position('none')

    a, b = np.polyfit(data["A"][data["x"]==label_r[counter]], data["F/A"][data["x"]==label_r[counter]],deg=1)
    y_est = a * data["A"][data["x"]==label_r[counter]] + b
    axs[2, column].plot(data["A"][data["x"]==label_r[counter]],
                        data["F/A"][data["x"]==label_r[counter]], 
                        "o", ms=2, color="tab:red")
    axs[2, column].plot(data["A"][data["x"]==label_r[counter]], y_est,
                        color="tab:blue", lw=1)
    axs[2, column].spines["right"].set_visible(False)
    axs[2, column].spines["top"].set_visible(False)
    axs[2, column].set_ylim([0, 10000000])
    axs[2, column].set_yticks(np.arange(0, 12000000, 10000000))
    axs[2, column].set_xlim([0, 1])
    axs[2, column].set_xticks(np.arange(0, 1.2, 1))
    #axs[2, column].yaxis.set_ticks_position('none')
    #axs[2, column].xaxis.set_ticks_position('none')
    counter = counter + 1

#axs[1,0].text(-1.5, -1.4, "Cm$^R$", ha="center")
#axs[1,0].text(0.2, -2.5, "pSC101 ori", ha="center")
#axs[1,1].text(-1.5, -1.6, "Cm$^R$", ha="center")
#axs[1,1].text(1.5, -2, "p15A ori", ha="center")
#axs[1,2].text(-1.2, -2.1, "Cm$^R$", ha="center")
#axs[1,2].text(1.5, -2, "pUC ori", ha="center")
axs[0,0].set_ylabel("A$_{600}$")
axs[1,0].set_ylabel("mRFP / A$_{600}$")
axs[2,0].set_ylabel("mRFP / A$_{600}$")
axs[0,0].set_xlabel("Hour")
axs[0,1].set_xlabel("Hour")
axs[0,2].set_xlabel("Hour")
axs[1,0].set_xlabel("Hour")
axs[1,1].set_xlabel("Hour")
axs[1,2].set_xlabel("Hour")
axs[2,0].set_xlabel("A$_{600}$")
axs[2,1].set_xlabel("A$_{600}$")
axs[2,2].set_xlabel("A$_{600}$")
#axbig.text(-0.11, 0.5, "(a)", transform=axbig.transAxes, fontsize=9, va='top', ha='right')
axs[0,0].text(-0.5, 0.5, "(a)", transform=axs[0,0].transAxes, fontsize=9, va='top', ha='right')
axs[1,0].text(-0.5, 0.5, "(b)", transform=axs[1,0].transAxes, fontsize=9, va='top', ha='right')
axs[2,0].text(-0.5, 0.5, "(c)", transform=axs[2,0].transAxes, fontsize=9, va='top', ha='right')
#axs[3,0].text(-0.5, 0.5, "(d)", transform=axs[3,0].transAxes, fontsize=9, va='top', ha='right')
#axs[4,0].text(-0.4, 0.5, "(e)", transform=axs[4,0].transAxes, fontsize=9, va='top', ha='right')

axs[row, 0].set_title(label="low copy", fontdict={"fontsize": 9})
axs[row, 1].set_title(label="medium copy", fontdict={"fontsize": 9})
axs[row, 2].set_title(label="high copy", fontdict={"fontsize": 9})

fig.tight_layout()
#fig.savefig("../figure/thesis-mrfp-graph.pdf")
fig.savefig("../figure/thesis-mrfp-graph.png", dpi=300, bbox_inches='tight', transparent=True)