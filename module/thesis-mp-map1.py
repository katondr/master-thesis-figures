#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:49:54 2021

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

fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, sharey="row",
                           facecolor='white', figsize=(6, 3))

counter = 0

features_2=[
    GraphicFeature(start=0, end=20, strand=+1, color="tab:grey"),
    GraphicFeature(start=77, end=117, strand=-1, color="tab:red"),
    GraphicFeature(start=141, end=195, strand=+1, color="tab:blue"),
    GraphicFeature(start=204, end=215, strand=+1, color="tab:blue"),
    GraphicFeature(start=222, end=944, strand=+1, color="tab:olive"),
    GraphicFeature(start=953, end=1032, strand=+1, color="tab:red"),
    GraphicFeature(start=1041, end=1081, strand=+1, color="tab:red"),
    GraphicFeature(start=1103, end=1160, strand=+1, color="tab:red"),
    GraphicFeature(start=1236, end=1255, strand=-1, color="tab:grey")
]

linear_features_2 = copy.deepcopy(features_2)

#linear_features_2[0].label = "VF2"
#linear_features_2[0].fontdict["fontsize"] = 9
#linear_features_2[2].label = "$\it{luxI}$ promoter"
#linear_features_2[2].fontdict["fontsize"] = 9
#linear_features_2[3].label = "BBa_B0034"
#linear_features_2[3].fontdict["fontsize"] = 9
linear_features_2[4].label = "EYFP"
linear_features_2[4].fontdict["fontsize"] = 9
#linear_features_2[5].label = "BBa_B0010"
#linear_features_2[5].fontdict["fontsize"] = 9
#linear_features_2[6].label = "BBa_B0012"
#linear_features_2[6].fontdict["fontsize"] = 9
#linear_features_2[8].label = "VR"
#linear_features_2[8].fontdict["fontsize"] = 9


plasmid_6_features = features_2.copy() + [
    GraphicFeature(start=1304, end=2122, strand=-1, color="tab:grey"),
    GraphicFeature(start=2480, end=3430, strand=-1, color="tab:purple")

]
plasmid_7_features = features_2.copy() + [
    GraphicFeature(start=1318, end=2133, strand=-1, color="tab:grey"),
    GraphicFeature(start=2465, end=3289, strand=-1, color="tab:purple")

]
plasmid_8_features = features_2.copy() + [
    GraphicFeature(start=1355, end=1943, strand=-1, color="tab:grey"),
    GraphicFeature(start=2241, end=2900, strand=-1, color="tab:green")

]

gs_2 = axs[0, 0].get_gridspec()
axs[0,0].remove()
axs[0,1].remove()
axs[0,2].remove()

record_2 = GraphicRecord(sequence_length=1255, features=linear_features_2)
axbig_2 = fig.add_subplot(gs_2[0,:])
record_2.plot(ax=axbig_2)
#axbig_2.set_xticks(np.arange(0, 1300, 300))
axbig_2.set_xticks([])

plasmid_6 = CircularGraphicRecord(sequence_length=4394, features=plasmid_6_features,
                                 top_position=1255/2)
plasmid_6.plot(ax=axs[1,0])

plasmid_7 = CircularGraphicRecord(sequence_length=3689, features=plasmid_7_features,
                                 top_position=1255/2)
plasmid_7.plot(ax=axs[1,1])

plasmid_8 = CircularGraphicRecord(sequence_length=3011, features=plasmid_8_features,
                                 top_position=1255/2)
plasmid_8.plot(ax=axs[1,2])

axs[1,0].text(-0.2, -2.4, "Km$^R$", ha="center")
axs[1,0].text(1.8, -1.3, "pSC101\nori", ha="center")

axs[1,1].text(-1, -2.2, "Km$^R$", ha="center")
axs[1,1].text(1.6, -1.6, "p15A\nori", ha="center")

axs[1,2].text(-1.2, -2, "Cm$^R$", ha="center")
axs[1,2].text(1.3, -2.2, "pUC\nori", ha="center")

axbig_2.text(-0.1, 0.5, "(a)", transform=axbig_2.transAxes, fontsize=9, va='top', ha='right')
axs[1,0].text(-0.4, 0.5, "(b)", transform=axs[1,0].transAxes, fontsize=9, va='top', ha='right')
#axs[2,0].text(-0.4, 0.5, "(c)", transform=axs[2,0].transAxes, fontsize=12, va='top', ha='right')
#axbig_2.text(-0.1, 0.5, "(d)", transform=axbig_2.transAxes, fontsize=12, va='top', ha='right')
#axs[4,0].text(-0.4, 0.5, "(e)", transform=axs[4,0].transAxes, fontsize=12, va='top', ha='right')

fig.tight_layout()
#fig.savefig("../figure/thesis-mp-map1.pdf")
fig.savefig("../figure/thesis-mp-map1.png", dpi=300, bbox_inches='tight', transparent=True)