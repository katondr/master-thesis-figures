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
features=[
    GraphicFeature(start=0, end=20, strand=+1, color="tab:grey"),
    GraphicFeature(start=77, end=117, strand=-1, color="tab:red"),
    GraphicFeature(start=274, end=1026, strand=-1, color="gold"),
    #GraphicFeature(start=1072, end=1105, strand=-1, color="tab:blue"),
    #GraphicFeature(start=1184, end=1228, strand=+1, color="tab:blue"),
    GraphicFeature(start=1245, end=1887, strand=+1, color="gold"),
    GraphicFeature(start=1896, end=1907, strand=+1, color="blue"),
    GraphicFeature(start=1914, end=2630, strand=+1, color="tab:olive"),
    GraphicFeature(start=2658, end=2715, strand=+1, color="tab:red"),
    GraphicFeature(start=2811, end=2830, strand=-1, color="tab:grey")
]
linear_features = copy.deepcopy(features)
#linear_features[0].label = "VF2"
#linear_features[0].fontdict["fontsize"] = 9
linear_features[2].label = "$\it{luxR}$"
linear_features[2].fontdict["fontsize"] = 9
#linear_features[3].label = "$\it{luxR}$ promoter"
#linear_features[3].fontdict["fontsize"] = 9
#linear_features[4].label = "$\it{luxI}$ promoter"
#linear_features[4].fontdict["fontsize"] = 9
linear_features[5-2].label = "$\it{luxI}$"
linear_features[5-2].fontdict["fontsize"] = 9
#linear_features[6].label = "BBa_B0034"
#linear_features[6].fontdict["fontsize"] = 9
linear_features[7-2].label = "sYFP2"
linear_features[7-2].fontdict["fontsize"] = 9
#linear_features[9].label = "VR"
#linear_features[9].fontdict["fontsize"] = 9
record = GraphicRecord(sequence_length=2830, features=linear_features)
plasmid_1_features = features.copy() + [
    GraphicFeature(start=2930, end=3518, strand=-1, color="tab:grey"),
    GraphicFeature(start=3816, end=4475, strand=-1, color="tab:green")
]
plasmid_2_features = features.copy() + [
    GraphicFeature(start=3277, end=4101, strand=+1, color="tab:grey"),
    GraphicFeature(start=4547, end=5209, strand=+1, color="tab:green")
]
plasmid_3_features = features.copy() + [
    GraphicFeature(start=3764, end=4714, strand=+1, color="tab:grey"),
    GraphicFeature(start=5030, end=5692, strand=+1, color="tab:green")
]
plasmid_1 = CircularGraphicRecord(sequence_length=4586, features=plasmid_1_features,
                                 top_position=2830/2)
plasmid_2 = CircularGraphicRecord(sequence_length=5288, features=plasmid_2_features,
                                 top_position=2830/2)
plasmid_3 = CircularGraphicRecord(sequence_length=5771, features=plasmid_3_features,
                                 top_position=2830/2)
#%%
#fig, axs = plt.subplots(111,
#                        sharex=False,
#                        sharey="row",
                        #facecolor='white',
                        #figsize=(4, 6))
#fig = plt.figure(figsize=(6, 2))
fig, _ = record.plot(figure_width=6, figure_height=1.25)
fig.set_xticks([])

fig.figure.tight_layout()
#fig.tight_layout(pad=0, w_pad=0.2, h_pad=0.2)
#fig.figure.savefig("../figure/thesis-luxri-34y-map.pdf")
fig.figure.savefig("../figure/thesis-luxri-34y-map.png", dpi=300, bbox_inches='tight', transparent=True)