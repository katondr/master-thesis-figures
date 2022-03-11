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
#%%
counter = 0
features=[
    GraphicFeature(start=0, end=20, strand=+1, color="tab:grey"),
    GraphicFeature(start=77, end=117, strand=-1, color="tab:red"),
    GraphicFeature(start=141, end=195, strand=+1, color="tab:blue"),
    GraphicFeature(start=204, end=215, strand=+1, color="tab:blue"),
    GraphicFeature(start=222, end=944, strand=+1, color="tab:olive"),
    GraphicFeature(start=953, end=1081, strand=+1, color="tab:red"),
    GraphicFeature(start=1223, end=1975, strand=-1, color="tab:orange"),
    GraphicFeature(start=2021, end=2054, strand=-1, color="tab:blue"),
    GraphicFeature(start=2123, end=2177, strand=+1, color="tab:blue"),
    GraphicFeature(start=2194, end=2836, strand=+1, color="tab:orange"),
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
#%%
record = GraphicRecord(sequence_length=3030, features=linear_features)
    
fig, _ = record.plot(figure_width=6, figure_height=1.5)
fig.set_xticks(np.arange(0, 3200, 750))

fig.figure.tight_layout()

#fig.savefig("../figure/thesis-34y-luxri.png", dpi=300, transparent=True)
#fig.figure.savefig("../figure/thesis-34y-luxri-map.pdf")
#fig.figure.savefig("../figure/thesis-34y-luxri-map.png", dpi=300, bbox_inches='tight', transparent=True)