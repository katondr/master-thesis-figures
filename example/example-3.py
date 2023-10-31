#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:51:08 2021

@author: Katon Dorojatun
"""
from setting import *
#%%
counter = 0
features=[
    GraphicFeature(start=0, end=20, strand=+1, color="tab:grey"),
    GraphicFeature(start=77, end=117, strand=-1, color="tab:red"),
    GraphicFeature(start=141, end=195, strand=+1, color="tab:blue"),
    GraphicFeature(start=204, end=215, strand=+1, color="tab:blue"),
    GraphicFeature(start=222, end=944, strand=+1, color="tab:cyan"),
    GraphicFeature(start=953, end=1081, strand=+1, color="tab:red"),
    GraphicFeature(start=1223, end=1975, strand=-1, color="tab:olive"),
    GraphicFeature(start=2021, end=2054, strand=-1, color="tab:blue"),
    GraphicFeature(start=2123, end=2177, strand=+1, color="tab:blue"),
    GraphicFeature(start=2194, end=2836, strand=+1, color="tab:olive"),
    GraphicFeature(start=2858, end=2915, strand=+1, color="tab:red"),
    GraphicFeature(start=3011, end=3030, strand=-1, color="tab:grey"),
]
linear_features = copy.deepcopy(features)

linear_features[4].label = "reporter"
linear_features[4].fontdict["fontsize"] = 9

linear_features[6].label = "gene A"
linear_features[6].fontdict["fontsize"] = 9

linear_features[9].label = "gene B"
linear_features[9].fontdict["fontsize"] = 9
#%%
record = GraphicRecord(sequence_length=3030, features=linear_features)
    
fig, _ = record.plot(figure_width=6, figure_height=1.5)
fig.set_xticks(np.arange(0, 3200, 750))

#fig.figure.tight_layout()
#fig.savefig("../figure/thesis-34y-luxri.png", dpi=300, transparent=True)
#fig.figure.savefig("../figure/thesis-34y-luxri-map.pdf")
fig.figure.savefig("figure/fig-3.png", dpi=300, bbox_inches='tight', transparent=False)
