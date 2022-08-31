#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:51:08 2021

@author: Katon Dorojatun
"""

from setting import *

#%%
data = pd.read_csv('../data/fl-210615.csv')
#time = data.t.unique().tolist()[:-1]
time_short = ["0", "", "4", "", "8", "", "12"]
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
#%%
fig, axs = plt.subplots(nrows=3, ncols=3, sharex="row", sharey="row",
                           facecolor='white', figsize=(5, 5.5))
width = 0.35
x = np.arange(len(time))
counter = 0
row = 0
for column in np.arange(len(label)):
    axs[row, column].errorbar(x, 
                            result["A_average"][result["x"] == label[counter]],
                            yerr=result["A_std"][result["x"] == label[counter]],
                            ecolor="tab:grey",  color="tab:orange",
                            capsize=2, capthick=0.5, lw=1, mec="tab:grey",
                            marker=".", markersize=1)
    axs[row, column].spines["right"].set_visible(False)
    axs[row, column].spines["top"].set_visible(False)
    axs[row, column].set_xticks(x)
    axs[row, column].set_xticklabels(time_short)
    axs[row, column].set_ylim([0, 1])
    axs[row, column].set_yticks(np.arange(0, 1.2, 1))
    #axs[row, column].yaxis.set_ticks_position('none')
    #axs[row, column].xaxis.set_ticks_position('none')
    
    axs[row+1, column].errorbar(x,
                            result["F/A_average"][result["x"] == label[counter]],
                            yerr= result["F/A_std"][result["x"] == label[counter]],
                            ecolor="tab:grey",  color="tab:olive",
                            capsize=2, capthick=0.5, lw=1, mec="tab:gray",
                            marker=".", markersize=1)
    axs[row+1, column].spines["right"].set_visible(False)
    axs[row+1, column].spines["top"].set_visible(False)
    axs[row+1, column].set_xticks(x)
    axs[row+1, column].set_xticklabels(time_short)
    axs[row+1, column].set_ylim([0, 120000000])
    axs[row+1, column].set_yticks(np.arange(0, 120000000, 100000000))
    #axs[row+1, column].yaxis.set_ticks_position('none')
    #axs[row+1, column].xaxis.set_ticks_position('none')

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
  
    axs[row+2, column].plot(x_axis, y_axis, "o", ms=2, color="tab:olive")
    axs[row+2, column].spines["right"].set_visible(False)
    axs[row+2, column].spines["top"].set_visible(False)
    axs[row+2, column].set_ylim([0, 120000000])
    axs[row+2, column].set_yticks(np.arange(0, 120000000, 100000000))
    axs[row+2, column].set_xlim([0, 1])
    axs[row+2, column].set_xticks(np.arange(0, 1.2, 1))
    #axs[row+2, column].yaxis.set_ticks_position('none')
    #axs[row+2, column].xaxis.set_ticks_position('none')
    def exp(x, a, b, y0):
        return (a * np.exp(b*x)) + y0
    #p0 = [70000, 9, min(y_axis)] # this is an mandatory initial guess
    popt, pcov = opt.curve_fit(exp, x_axis, y_axis, method="dogbox")
    x_fit = np.linspace(min(x_axis), max(x_axis), num= 300)
    y_fit = exp(x_fit, *popt)
    axs[row+2, column].plot(x_fit, y_fit, color="tab:blue", lw=1)
    print(popt)
    counter = counter + 1
    
axs[row,0].set_ylabel("A$_{600}$")
axs[row+1,0].set_ylabel("EYFP / A$_{600}$")
axs[row+2,0].set_ylabel("EYFP / A$_{600}$")
axs[row,0].set_xlabel("Hour")
axs[row,1].set_xlabel("Hour")
axs[row,2].set_xlabel("Hour")
axs[row+1,0].set_xlabel("Hour")
axs[row+1,1].set_xlabel("Hour")
axs[row+1,2].set_xlabel("Hour")
axs[row+2,0].set_xlabel("A$_{600}$")
axs[row+2,1].set_xlabel("A$_{600}$")
axs[row+2,2].set_xlabel("A$_{600}$")

axs[row,0].text(-0.5, 0.5, "(a)", transform=axs[row,0].transAxes, fontsize=9, va='top', ha='right')
axs[row+1,0].text(-0.5, 0.5, "(b)", transform=axs[row+1,0].transAxes, fontsize=9, va='top', ha='right')
axs[row+2,0].text(-0.5, 0.5, "(c)", transform=axs[row+2,0].transAxes, fontsize=9, va='top', ha='right')
#axs[row+1,0].text(-0.4, 0.5, "(d)", transform=axs[3,0].transAxes, fontsize=9, va='top', ha='right')
#axs[row+2,0].text(-0.4, 0.5, "(e)", transform=axs[4,0].transAxes, fontsize=9, va='top', ha='right')

axs[row, 0].set_title(label="low copy", fontdict={"fontsize": 9})
axs[row, 1].set_title(label="medium copy", fontdict={"fontsize": 9})
axs[row, 2].set_title(label="high copy", fontdict={"fontsize": 9})

fig.tight_layout()
#fig.savefig("../figure/thesis-34y-luxri.png", dpi=300, transparent=True)
fig.savefig("../figure/thesis-34y-luxri-graph.pdf")
fig.savefig("../figure/thesis-34y-luxri-graph.png", dpi=300, bbox_inches='tight', transparent=True)