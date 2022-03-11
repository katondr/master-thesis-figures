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
#linear_features[0].label = "VF2"
#linear_features[0].fontdict["fontsize"] = 9
#linear_features[3].label = "B0034"
#linear_features[3].fontdict["fontsize"] = 9
#linear_features[2].label = "lac promoter"
#linear_features[2].fontdict["fontsize"] = 9
linear_features[4].label = "mRFP"
linear_features[4].fontdict["fontsize"] = 9
#linear_features[5].label = "nrB T1"
linear_features[5].fontdict["fontsize"] = 9
#linear_features[7].label = "VR"
#linear_features[7].fontdict["fontsize"] = 9
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
fig, _ = record.plot(figure_width=6, figure_height=1.25)
fig.set_xticks(np.arange(0, 1383, 340))
fig.figure.tight_layout()
#fig.figure.savefig("../figure/thesis-mrfp-map.pdf")
fig.figure.savefig("../figure/thesis-mrfp-map.png", dpi=300, bbox_inches='tight', transparent=True)