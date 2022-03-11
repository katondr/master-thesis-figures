#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:46:04 2021

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

data = pd.read_csv('../data/fl-210608.csv')
time = ["0", "2", "4", "6", "8", "10", "12"]
time_short = ["0", "", "4", "", "8", "", "12"]
label = data.x.unique().tolist()
abm = data["A"][data["t"] == "blc"].mean()
bfbm = data["BF"][data["t"] == "blc"].mean()
yfbm = data["YF"][data["t"] == "blc"].mean()
data["A"] -= abm
data["BF"] -= bfbm
data["YF"] -= yfbm
data["BF/A"] = data["BF"]/data["A"]
data["YF/A"] = data["YF"]/data["A"]
data = data.loc[data["t"] != "blc"]
result = pd.DataFrame({"x":[], "t":[], "A_average":[], "A_std":[],
                       "BF/A_average":[], "YF/A_average":[], "BF/A_std":[], "YF/A_std":[]})

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
    result.loc[index, "BF/A_average"] = data["BF/A"][
        (data["x"]==row["x"])&(data["t"]==row["t"])
    ].mean()
for index, row in result.iterrows():
    result.loc[index, "YF/A_average"] = data["YF/A"][
        (data["x"]==row["x"])&(data["t"]==row["t"])
    ].mean()
for index, row in result.iterrows():
    result.loc[index, "BF/A_std"] = data["BF/A"][
        (data["x"]==row["x"])&(data["t"]==row["t"])
    ].std()
for index, row in result.iterrows():
    result.loc[index, "YF/A_std"] = data["YF/A"][
        (data["x"]==row["x"])&(data["t"]==row["t"])
    ].std()

label_1 = label
    
fig, axs = plt.subplots(nrows=5, ncols=4, sharex=False, sharey="row",
                           facecolor='white', figsize=(5.5, 8))

x = np.arange(len(time))

counter = 0
for column in np.arange(len(label_1)):
    axs[0, column].errorbar(x,
                            result["A_average"][result["x"] == label_1[counter]],
                            yerr=result["A_std"][result["x"] == label_1[counter]],
                            ecolor="tab:grey",  color="tab:orange",
                            capsize=2, capthick=0.5, lw=1, mec="tab:grey",
                            marker=".", markersize=1)
    axs[0, column].spines["right"].set_visible(False)
    axs[0, column].spines["top"].set_visible(False)
    axs[0, column].set_xticks(x)
    axs[0, column].set_xticklabels(time_short)
    axs[0, column].set_ylim([0, 1])
    axs[0, column].set_yticks(np.arange(0, 1.2, 1))
    #axs[0, column].set_yticklabels([0, 0.3, 0.6])
    #axs[0, column].yaxis.set_ticks_position('none')
    #axs[0, column].xaxis.set_ticks_position('none')
    
    axs[1, column].errorbar(x,
                            result["BF/A_average"][result["x"] == label_1[counter]],
                            yerr= result["BF/A_std"][result["x"] == label_1[counter]],
                            ecolor="tab:grey",  color="tab:blue",
                            capsize=2, capthick=0.5, lw=1, mec="tab:gray",
                            marker=".", markersize=1)
    axs[1, column].spines["right"].set_visible(False)
    axs[1, column].spines["top"].set_visible(False)
    axs[1, column].set_xticks(x)
    axs[1, column].set_xticklabels(time_short)
    axs[1, column].set_ylim([-100000000, 300000000])
    axs[1, column].set_yticks(np.arange(0, 210000000, 200000000))
    #axs[1, column].yaxis.set_ticks_position('none')
    #axs[1, column].xaxis.set_ticks_position('none')
        
    selector = data["x"]==label_1[counter]   
    temp_data = result["A_average"][result["x"]==label_1[counter]]
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
        item_masking = (data["t"] == time[item]) & (selector)
        bool_masking = bool_masking | item_masking    
            
    x_axis = data["A"][(selector) & bool_masking]
    y_axis = data["BF/A"][(selector) & bool_masking]
    
    axs[2, column].plot(data["A"][selector],
                        data["BF/A"][selector],
                        "o", ms=2, color="tab:blue")
    axs[2, column].spines["right"].set_visible(False)
    axs[2, column].spines["top"].set_visible(False)
    axs[2, column].set_ylim([-100000000, 300000000])
    axs[2, column].set_yticks(np.arange(0, 210000000, 200000000))
    axs[2, column].set_xlim([0, 1])
    axs[2, column].set_xticks(np.arange(0, 1.2, 1))
    #axs[2, column].set_xticklabels([0, 0.4, 0.8])
    #axs[2, column].yaxis.set_ticks_position('none')
    #axs[2, column].xaxis.set_ticks_position('none')
    
    def exp(x, a, b, y0):
        return a * np.exp(b*x) + y0
    p0 = [8000000, 55, min(y_axis)] # this is an mandatory initial guess
    popt, pcov = opt.curve_fit(exp, x_axis, y_axis, p0, method="dogbox")
    x_fit = np.linspace(min(x_axis), max(x_axis), num= 300)
    y_fit = exp(x_fit, *popt)       
    axs[2, column].plot(x_fit, y_fit, color="tab:red")

    y_axis = data["YF/A"][(selector) & bool_masking]
    
    axs[3, column].errorbar(x,
                            result["YF/A_average"][result["x"] == label_1[counter]],
                            yerr= result["YF/A_std"][result["x"] == label_1[counter]],
                            ecolor="tab:grey",  color="tab:olive",
                            capsize=2, capthick=0.5, lw=1, mec="tab:gray",
                            marker=".", markersize=1)
    axs[3, column].spines["right"].set_visible(False)
    axs[3, column].spines["top"].set_visible(False)
    axs[3, column].set_xticks(x)
    axs[3, column].set_xticklabels(time_short)
    axs[3, column].set_ylim([-5000000, 45000000])
    axs[3, column].set_yticks(np.arange(0, 41000000, 40000000))
    #axs[3, column].yaxis.set_ticks_position('none')
    #axs[3, column].xaxis.set_ticks_position('none')

    axs[4, column].plot(data["A"][data["x"]==label_1[counter]],
                        data["YF/A"][data["x"]==label_1[counter]],
                        "o", ms=2, color="tab:olive")
    axs[4, column].spines["right"].set_visible(False)
    axs[4, column].spines["top"].set_visible(False)
    axs[4, column].set_ylim([-5000000, 45000000])
    axs[4, column].set_yticks(np.arange(0, 41000000, 40000000))
    axs[4, column].set_xlim([0, 1])
    axs[4, column].set_xticks(np.arange(0, 1.2, 1))
    #axs[4, column].yaxis.set_ticks_position('none')
    #axs[4, column].xaxis.set_ticks_position('none')
    
    def exp(x, a, b, y0):
        return a * np.exp(b*x) + y0
    #p0 = [70000, 9, min(y_axis)] # this is an mandatory initial guess
    popt, pcov = opt.curve_fit(exp, x_axis, y_axis, method="dogbox")
    x_fit = np.linspace(min(x_axis), max(x_axis), num= 300)
    y_fit = exp(x_fit, *popt)   
    
    axs[4, column].plot(x_fit, y_fit, color="tab:blue")
    counter = counter + 1

axs[0,0].set_ylabel("A$_{600}$")
axs[1,0].set_ylabel("tagBFP / A$_{600}$")
axs[2,0].set_ylabel("tagBFP / A$_{600}$")
axs[3,0].set_ylabel("EYFP / A$_{600}$")
axs[4,0].set_ylabel("EYFP / A$_{600}$")
axs[0,0].set_xlabel("Hour")
axs[0,1].set_xlabel("Hour")
axs[0,2].set_xlabel("Hour")
axs[0,3].set_xlabel("Hour")
axs[1,0].set_xlabel("Hour")
axs[1,1].set_xlabel("Hour")
axs[1,2].set_xlabel("Hour")
axs[1,3].set_xlabel("Hour")
axs[2,0].set_xlabel("A$_{600}$")
axs[2,1].set_xlabel("A$_{600}$")
axs[2,2].set_xlabel("A$_{600}$")
axs[2,3].set_xlabel("A$_{600}$")
axs[3,0].set_xlabel("Hour")
axs[3,1].set_xlabel("Hour")
axs[3,2].set_xlabel("Hour")
axs[3,3].set_xlabel("Hour")
axs[4,0].set_xlabel("A$_{600}$")
axs[4,1].set_xlabel("A$_{600}$")
axs[4,2].set_xlabel("A$_{600}$")
axs[4,3].set_xlabel("A$_{600}$")

fig_label = ["(a)", "(b)", "(c)", "(d)", "(e)"]
for row in np.arange(len(axs)):
    axs[row,0].text(-0.5, 0.5, fig_label[row], transform=axs[row,0].transAxes, fontsize=9, va='top', ha='right')

axs[0,0].set_title(label="4C-luxRI-34BFP\n3K-P$_{luxI}$-34EYFPt", fontdict={"fontsize": 9})    
axs[0,1].set_title(label="3C-luxRI-34BFP\n4K-P$_{luxI}$-34EYFPt", fontdict={"fontsize": 9}) 
axs[0,2].set_title(label="1C-luxRI-34BFP\n4K-P$_{luxI}$-34EYFPt", fontdict={"fontsize": 9}) 
axs[0,3].set_title(label="1C-luxRI-34BFP\n3K-P$_{luxI}$-34EYFPt", fontdict={"fontsize": 9}) 
    
fig.tight_layout()
#fig.savefig("../figure/thesis-mp2.pdf")
fig.savefig("../figure/thesis-mp2.png", dpi=300, bbox_inches='tight', transparent=True)