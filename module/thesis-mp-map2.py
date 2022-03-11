#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:50:19 2021

@author: tony
"""
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

fig, axs = plt.subplots(nrows=3, ncols=3, sharex=False, sharey="row",
                           facecolor='white', figsize=(6, 4))

counter = 0
features_1=[
    GraphicFeature(start=0, end=20, strand=+1, color="tab:grey"),
    GraphicFeature(start=77, end=117, strand=-1, color="tab:red"),
    GraphicFeature(start=274, end=1026, strand=-1, color="gold"),
    GraphicFeature(start=1072, end=1105, strand=-1, color="tab:blue"),
    GraphicFeature(start=1184, end=1228, strand=+1, color="tab:blue"),
    GraphicFeature(start=1245, end=1887, strand=+1, color="gold"),
    GraphicFeature(start=1909, end=1966, strand=+1, color="tab:red"),
    GraphicFeature(start=2042, end=2061, strand=-1, color="tab:grey")
]

linear_features_1 = copy.deepcopy(features_1)
#linear_features_1[0].label = "VF2"
#linear_features_1[0].fontdict["fontsize"] = 9
linear_features_1[2].label = "$\it{luxR}$"
linear_features_1[2].fontdict["fontsize"] = 9
#linear_features_1[3].label = "$\it{luxR}$ promoter"
#linear_features_1[3].fontdict["fontsize"] = 9
#linear_features_1[4].label = "$\it{luxI}$ promoter"
#linear_features_1[4].fontdict["fontsize"] = 9
linear_features_1[5].label = "$\it{luxI}$"
linear_features_1[5].fontdict["fontsize"] = 9
#linear_features_1[7].label = "VR"
#linear_features_1[7].fontdict["fontsize"] = 9

plasmid_1_features = features_1.copy() + [
    GraphicFeature(start=2161, end=2749, strand=-1, color="tab:grey"),
    GraphicFeature(start=3047, end=3706, strand=-1, color="tab:green")

]
plasmid_2_features = features_1.copy() + [
    GraphicFeature(start=2508, end=3332, strand=+1, color="tab:grey"),
    GraphicFeature(start=3778, end=4440, strand=+1, color="tab:green")

]
plasmid_3_features = features_1.copy() + [
    GraphicFeature(start=2995, end=3945, strand=+1, color="tab:grey"),
    GraphicFeature(start=4261, end=4923, strand=+1, color="tab:green")

]
plasmid_4_features = features_1.copy() + [
    GraphicFeature(start=2110, end=2928, strand=-1, color="tab:grey"),
    GraphicFeature(start=3286, end=4236, strand=-1, color="tab:purple")

]
plasmid_5_features = features_1.copy() + [
    GraphicFeature(start=2124, end=2939, strand=-1, color="tab:grey"),
    GraphicFeature(start=3271, end=4095, strand=-1, color="tab:purple")

]

gs_1 = axs[0, 0].get_gridspec()
axs[0,0].remove()
axs[0,1].remove()
axs[0,2].remove()

record_1 = GraphicRecord(sequence_length=2061, features=linear_features_1)
axbig_1 = fig.add_subplot(gs_1[0,:])
record_1.plot(ax=axbig_1)
#axbig_1.set_xticks(np.arange(0, 2100, 500))
axbig_1.set_xticks([])

plasmid_1 = CircularGraphicRecord(sequence_length=3817, features=plasmid_1_features,
                                 top_position=2061/2)
plasmid_1.plot(ax=axs[1,2])

plasmid_2 = CircularGraphicRecord(sequence_length=4519, features=plasmid_2_features,
                                 top_position=2061/2)
plasmid_2.plot(ax=axs[1,1])

plasmid_3 = CircularGraphicRecord(sequence_length=5002, features=plasmid_3_features,
                                 top_position=2061/2)
plasmid_3.plot(ax=axs[1,0])

plasmid_4 = CircularGraphicRecord(sequence_length=5200, features=plasmid_4_features,
                                 top_position=2061/2)
plasmid_4.plot(ax=axs[2,0])

plasmid_5 = CircularGraphicRecord(sequence_length=4495, features=plasmid_5_features,
                                 top_position=2061/2)
plasmid_5.plot(ax=axs[2,1])

axs[2,2].remove()

axs[1,0].text(-1.5, -1.6, "Cm$^R$", ha="center")
axs[1,0].text(0.2, -2.5, "pSC101 ori", ha="center")

axs[1,1].text(-1.5, -1.8, "Cm$^R$", ha="center")
axs[1,1].text(0.8, -2.4, "p15A ori", ha="center")

axs[1,2].text(-1.2, -2.2, "Cm$^R$", ha="center")
axs[1,2].text(1.4, -2.2, "pUC ori", ha="center")

axs[2,0].text(-0.2, -2.5, "Km$^R$", ha="center")
axs[2,0].text(1.9, -1.7, "pSC101\nori", ha="center")

axs[2,1].text(-0.8, -2.5, "Km$^R$", ha="center")
axs[2,1].text(1.6, -2, "p15A\nori", ha="center")

axbig_1.text(-0.1, 0.5, "(a)", transform=axbig_1.transAxes, fontsize=9, va='top', ha='right')
axs[1,0].text(-0.4, 0.5, "(b)", transform=axs[1,0].transAxes, fontsize=9, va='top', ha='right')
axs[2,0].text(-0.4, 0.5, "(c)", transform=axs[2,0].transAxes, fontsize=9, va='top', ha='right')
#axbig_2.text(-0.1, 0.5, "(d)", transform=axbig_2.transAxes, fontsize=12, va='top', ha='right')
#axs[4,0].text(-0.4, 0.5, "(e)", transform=axs[4,0].transAxes, fontsize=12, va='top', ha='right')

fig.tight_layout()
#fig.savefig("../figure/thesis-mp-map2.pdf")
fig.savefig("../figure/thesis-mp-map2.png", dpi=300, bbox_inches='tight', transparent=True)