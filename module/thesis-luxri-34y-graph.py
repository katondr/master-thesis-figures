#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 11:46:31 2021

@author: tony
"""
#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as opt
import copy
from matplotlib.ticker import ScalarFormatter
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.sans-serif"] = "Times New Roman"
plt.rcParams['font.size'] = "9"
#%%
data = pd.read_csv('../data/fl-210526.csv')
time = ["0", "2", "4", "6", "8", "10", "12"]
label = data.x.unique().tolist()
abm = data["A"][data["t"] == "blc"].mean()
fbm = data["F"][data["t"] == "blc"].mean()
data["A"] -= abm
data["F"] -= fbm
data["F/A"] = data["F"]/data["A"]
data = data.loc[data["t"] != "blc"]
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
#%%
fig, axs = plt.subplots(nrows=3,
                        ncols=3,
                        sharex=False,
                        sharey="row",
                        facecolor='white',
                        figsize=(5, 5.5))
width = 0.35
x = np.arange(len(time))
counter = 0
row=0
for column in np.arange(len(label)):
    strain_selector = result["x"] == label[counter]
    axs[0, column].errorbar(x,
                            result["A_average"][strain_selector],
                            yerr=result["A_std"][strain_selector],
                            ecolor="tab:grey",  color="tab:orange",
                            capsize=2, capthick=0.5, lw=1, mec="tab:grey",
                            marker=".", markersize=1)
    axs[0, column].spines["right"].set_visible(False)
    axs[0, column].spines["top"].set_visible(False)
    axs[0, column].set_xticks(x)
    axs[0, column].set_xticklabels(time)
    axs[0, column].set_ylim([0, 1])
    axs[0, column].set_yticks(np.arange(0, 1.2, 0.5))
    #axs[1, column].yaxis.set_ticks_position('none')
    #axs[1, column].xaxis.set_ticks_position('none')

    axs[1, column].errorbar(x,
                            result["F/A_average"][strain_selector],
                            yerr= result["F/A_std"][strain_selector],
                            ecolor="tab:grey",  color="tab:olive",
                            capsize=2, capthick=0.5, lw=1, mec="tab:gray",
                            marker=".", markersize=1)
    axs[1, column].spines["right"].set_visible(False)
    axs[1, column].spines["top"].set_visible(False)
    axs[1, column].set_xticks(x)
    axs[1, column].set_xticklabels(time)
    axs[1, column].set_yticks(np.arange(0, 120000000, 50000000))
    #axs[2, column].yaxis.set_ticks_position('none')
    #axs[2, column].xaxis.set_ticks_position('none')

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
    axs[2, column].plot(x_axis, y_axis, "o", ms=2, color="tab:olive")
    axs[2, column].spines["right"].set_visible(False)
    axs[2, column].spines["top"].set_visible(False)
    axs[2, column].set_yticks(np.arange(0, 120000000, 50000000))
    axs[2, column].set_xlim([0, 1])
    axs[2, column].set_xticks(np.arange(0, 1.2, 0.5))
    #axs[3, column].yaxis.set_ticks_position('none')
    #axs[3, column].xaxis.set_ticks_position('none')
    def exp(x, a, b, y0):
        return a * np.exp(b*x) + y0
    #p0 = [70000, 9, min(y_axis)] # this is an mandatory initial guess
    popt, pcov = opt.curve_fit(exp, x_axis, y_axis, method="dogbox")
    x_fit = np.linspace(min(x_axis), max(x_axis), num= 300)
    y_fit = exp(x_fit, *popt)
    axs[2, column].plot(x_fit, y_fit, color="tab:blue", lw=1)
    counter = counter + 1
#axbig.set_xticks(np.arange(0, 2800, 900))
#axs[1,0].text(-1.5, -1.8, "Cm$^R$", ha="center")
#axs[1,0].text(0.2, -2.5, "pSC101 ori", ha="center")
#axs[1,1].text(-1.4, -2, "Cm$^R$", ha="center")
#axs[1,1].text(1, -2.4, "p15A ori", ha="center")
#axs[1,2].text(-1, -2.4, "Cm$^R$", ha="center")
#axs[1,2].text(1, -2.4, "pUC ori", ha="center")
axs[0,0].set_ylabel("A$_{600}$")
axs[1,0].set_ylabel("sYFP2 / A$_{600}$")
axs[2,0].set_ylabel("sYFP2 / A$_{600}$")
axs[0,0].set_xlabel("Hour")
axs[0,1].set_xlabel("Hour")
axs[0,2].set_xlabel("Hour")
axs[1,0].set_xlabel("Hour")
axs[1,1].set_xlabel("Hour")
axs[1,2].set_xlabel("Hour")
axs[2,0].set_xlabel("A$_{600}$")
axs[2,1].set_xlabel("A$_{600}$")
axs[2,2].set_xlabel("A$_{600}$")
axs[0,0].text(-0.6, 0.5, "(a)", transform=axs[0,0].transAxes, fontsize=9, va='top', ha='right')
axs[1,0].text(-0.6, 0.5, "(b)", transform=axs[1,0].transAxes, fontsize=9, va='top', ha='right')
axs[2,0].text(-0.6, 0.5, "(c)", transform=axs[2,0].transAxes, fontsize=9, va='top', ha='right')
axs[row, 0].set_title(label="low copy", fontdict={"fontsize": 9})
axs[row, 1].set_title(label="medium copy", fontdict={"fontsize": 9})
axs[row, 2].set_title(label="high copy", fontdict={"fontsize": 9})
fig.tight_layout()
#fig.tight_layout(pad=0, w_pad=0.2, h_pad=0.2)
#fig.savefig("../figure/thesis-luxri-34y-graph.pdf")
fig.savefig("../figure/thesis-luxri-34y-graph.png", dpi=300, bbox_inches='tight', transparent=True)