import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy
from matplotlib.ticker import ScalarFormatter
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.sans-serif"] = "Times New Roman"
plt.rcParams['font.size'] = "9"

seq_file = open('data/fig-7-seq.txt')
sequence = seq_file.read()
seq_file.close()

features=[
    GraphicFeature(start=1, end=264, strand=+1, color="tab:olive"),
    GraphicFeature(start=397, end=747, strand=+1, color="tab:olive"),
    GraphicFeature(start=929, end=2158, strand=+1, color="tab:olive"),
    GraphicFeature(start=2201, end=4315, strand=+1, color="gold"),
    GraphicFeature(start=4261, end=4971, strand=+1, color="tab:olive"),
    GraphicFeature(start=4968, end=6347, strand=+1, color="tab:olive"),
    GraphicFeature(start=6461, end=6856, strand=+1, color="tab:olive"),
    GraphicFeature(start=7169, end=7515, strand=+1, color="tab:olive")
]
features[0].label = "slr0374"
features[0].fontdict["fontsize"] = 9
features[1].label = "slr0376"
features[1].fontdict["fontsize"] = 9
features[2].label = "slr0377"
features[2].fontdict["fontsize"] = 9
features[3].label = "slr0378"
features[3].fontdict["fontsize"] = 9
features[4].label = "slr0379"
features[4].fontdict["fontsize"] = 9
features[5].label = "slr0380"
features[5].fontdict["fontsize"] = 9
features[6].label = "slr0381"
features[6].fontdict["fontsize"] = 9
features[7].label = "fabG"
features[7].fontdict["fontsize"] = 9
gene = GraphicRecord(sequence=sequence, sequence_length=7515, features=features)
genes = gene.crop((1327, 5715))
fig, _ = genes.plot(figure_width=6, figure_height=1.5)
#fig.set_xticks(np.arange(0, 2800, 900))

#fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(5, 3))
#plasmid_a.plot(ax=ax1)
#plasmid_b.plot(ax=ax2)
#ax1.text(0.5, 0.5, "pSVAK\n5617 bp", transform=ax1.transAxes, fontsize=9, va='top', ha='center')
#ax1.text(0.75, 0.8, "mRFP", transform=ax1.transAxes, fontsize=9, va='top', ha='left')
#ax1.text(0.87, 0.35, "pUC ori", transform=ax1.transAxes, fontsize=9, va='top', ha='left')
#ax1.text(0.55, 0.09, "Km$^R$", transform=ax1.transAxes, fontsize=9, va='top', ha='left')
#ax1.text(0.18, 0.25, "ORF1", transform=ax1.transAxes, fontsize=9, va='top', ha='right')
#ax1.text(0.12, 0.65, "ORF2", transform=ax1.transAxes, fontsize=9, va='top', ha='right')
#ax1.text(0, 1, "(a)", transform=ax1.transAxes, fontsize=9, va='top', ha='right')
#ax2.text(0.5, 0.5, "pSVCK\n6491 bp", transform=ax2.transAxes, fontsize=9, va='top', ha='center')
#ax2.text(0.7, 0.82, "mRFP", transform=ax2.transAxes, fontsize=9, va='top', ha='left')
#ax2.text(0.9, 0.45, "pUC ori", transform=ax2.transAxes, fontsize=9, va='top', ha='left')
#ax2.text(0.7, 0.15, "Km$^R$", transform=ax2.transAxes, fontsize=9, va='top', ha='left')
#ax2.text(0.09, 0.52, "ORFb", transform=ax2.transAxes, fontsize=9, va='top', ha='right')
#ax2.text(0, 1, "(b)", transform=ax2.transAxes, fontsize=9, va='top', ha='right')
fig.set_xticks(np.arange(1327, 6000, 1000))
fig.set_xticklabels(np.arange(0, 5000, 1000))
fig.figure.tight_layout()
#fig.figure.savefig("../figure/thesis-slr0378-map.pdf")
fig.figure.savefig("figure/fig-7.png", dpi=300, bbox_inches='tight', transparent=True)
