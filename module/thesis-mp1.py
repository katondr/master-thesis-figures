#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 11:53:05 2021

@author: tony
"""
#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as opt
import copy
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.sans-serif"] = "Times New Roman"
plt.rcParams['font.size'] = "9"
#%%
data = pd.read_csv('../data/fl-210411.csv')
constant = pd.read_csv('../data/constant.csv')
time = ["0", "2", "4", "6", "8", "10", "12"]
time_short = ["0", "", "4", "", "8", "", "12"]
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
    
label_1 = label[3:]
#%%
fig, axs = plt.subplots(nrows=6, ncols=3, sharex="col", sharey="col",
                           facecolor='white', figsize=(5.5, 8))

x = np.arange(len(time))
width = 0.4
counter = 0
#row = 0
column = 0
for row in np.arange(len(label_1)):
    axs[row, column].errorbar(x,
                            result["A_average"][result["x"] == label_1[counter]],
                            yerr=result["A_std"][result["x"] == label_1[counter]],
                            ecolor="tab:grey",  color="tab:orange",
                            capsize=2, capthick=0.5, lw=1, mec="tab:grey",
                            marker=".", markersize=1)
    axs[row, column].spines["right"].set_visible(False)
    axs[row, column].spines["top"].set_visible(False)
    axs[row, column].set_xticks(x)
    axs[row, column].set_xticklabels(time_short)
    #axs[0, column].set_ylim([0, 1])
    axs[row, column].set_ylim([0, 0.6])
    #axs[0, column].set_yticks(np.arange(0, 1.2, 1))
    axs[row, column].set_yticks(np.arange(0, 0.7, 0.6))
    #axs[row, column].yaxis.set_ticks_position('none')
    #axs[row, column].xaxis.set_ticks_position('none')
    axs[row, column].set_ylabel("A$_{600}$")
    
    axs[row, column+1].errorbar(x,
                            result["F/A_average"][result["x"] == label_1[counter]],
                            yerr= result["F/A_std"][result["x"] == label_1[counter]],
                            ecolor="tab:grey",  color="tab:olive",
                            capsize=2, capthick=0.5, lw=1, mec="tab:gray",
                            marker=".", markersize=1)
    axs[row, column+1].spines["right"].set_visible(False)
    axs[row, column+1].spines["top"].set_visible(False)
    axs[row, column+1].set_xticks(x)
    axs[row, column+1].set_xticklabels(time_short)
    #axs[1, column].set_ylim([0, 100000000])
    axs[row, column+1].set_ylim([0, 50000000])
    #axs[1, column].set_yticks(np.arange(0, 110000000, 100000000))
    axs[row, column+1].set_yticks(np.arange(0, 60000000, 50000000))
    #axs[row, column+1].yaxis.set_ticks_position('none')
    #axs[row, column+1].xaxis.set_ticks_position('none')
    axs[row, column+1].set_ylabel("EYFP / A$_{600}$")
    
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
    y_axis = data["F/A"][(selector) & bool_masking]
    axs[row, column+2].plot(x_axis, y_axis, "o", ms=2, color="tab:olive")
    axs[row, column+2].spines["right"].set_visible(False)
    axs[row, column+2].spines["top"].set_visible(False)
    #axs[2, column].set_ylim([0, 100000000])
    axs[row, column+2].set_ylim([0, 50000000])
    #axs[2, column].set_yticks(np.arange(0, 110000000, 100000000))
    axs[row, column+2].set_yticks(np.arange(0, 60000000, 50000000))
    #axs[2, column].set_xlim([0, 1])
    axs[row, column+2].set_xlim([0, 0.6])
    #axs[2, column].set_xticks(np.arange(0, 1.2, 1))
    axs[row, column+2].set_xticks(np.arange(0, 0.7, 0.6))
    #axs[row, column+2].yaxis.set_ticks_position('none')
    #axs[row, column+2].xaxis.set_ticks_position('none')
    axs[row, column+2].set_ylabel("EYFP / A$_{600}$")
    
    def exp(x, a, b, y0):
        return a * np.exp(b*x) + y0
    #p0 = [70000, 9, min(y_axis)] # this is an mandatory initial guess
    popt, pcov = opt.curve_fit(exp, x_axis, y_axis, method="dogbox")
    x_fit = np.linspace(min(x_axis), max(x_axis), num= 300)
    y_fit = exp(x_fit, *popt)
    axs[row, column+2].plot(x_fit, y_fit, color="tab:blue", lw=1)
    counter = counter + 1
"""

axs[1,0].set_ylabel("EYFP / A$_{600}$")
axs[2,0].set_ylabel("EYFP / A$_{600}$")
axs[0,0].set_xlabel("Hour")
axs[0,1].set_xlabel("Hour")
axs[0,2].set_xlabel("Hour")
axs[0,3].set_xlabel("Hour")
axs[0,4].set_xlabel("Hour")
axs[0,5].set_xlabel("Hour")
axs[1,0].set_xlabel("Hour")
axs[1,1].set_xlabel("Hour")
axs[1,2].set_xlabel("Hour")
axs[1,3].set_xlabel("Hour")
axs[1,4].set_xlabel("Hour")
axs[1,5].set_xlabel("Hour")
axs[2,0].set_xlabel("A$_{600}$")
axs[2,1].set_xlabel("A$_{600}$")
axs[2,2].set_xlabel("A$_{600}$")
axs[2,3].set_xlabel("A$_{600}$")
axs[2,4].set_xlabel("A$_{600}$")
axs[2,5].set_xlabel("A$_{600}$")
"""
fig_label = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
for row in np.arange(len(axs)):
    axs[row,0].text(-0.5, 0.5, fig_label[row], transform=axs[row,0].transAxes, fontsize=9, va='top', ha='right')

#axs[0,0].set_title(label="4C-luxRI\n3K-P$_{luxI}$-34Yt", fontdict={"fontsize": 9})    
#axs[0,1].set_title(label="4K-luxRI\n1C-P$_{luxI}$-34Yt", fontdict={"fontsize": 9}) 
#axs[0,2].set_title(label="3C-luxRI\n4K-P$_{luxI}$-34Yt", fontdict={"fontsize": 9}) 
#axs[0,3].set_title(label="3K-luxRI\n1C-P$_{luxI}$-34Yt", fontdict={"fontsize": 9}) 
#axs[0,4].set_title(label="1C-luxRI\n4K-P$_{luxI}$-34Yt", fontdict={"fontsize": 9}) 
#axs[0,5].set_title(label="1C-luxRI\n3K-P$_{luxI}$-34Yt", fontdict={"fontsize": 9}) 
#axs[0,2].text(1, 0.5, "low copy\nquorum\nsensing\n\nmed copy\nreporter", transform=axs[0,2].transAxes, fontsize=9, va='center', ha='left')
#axs[0,2].text(1, 0.5, "low copy\nquorum\nsensing\n\nmed copy\nreporter", transform=axs[0,2].transAxes, fontsize=9, va='center', ha='left')

#axs[0,0].get_shared_x_axes().
fig.tight_layout()
#fig.savefig("../figure/thesis-mp1.pdf")
fig.savefig("../figure/thesis-mp1.png", dpi=300, bbox_inches='tight', transparent=True)