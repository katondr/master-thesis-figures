import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy
from matplotlib.ticker import ScalarFormatter
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.sans-serif"] = "Times New Roman"
plt.rcParams['font.size'] = "9"

file_a = open('fig-6-seq-1.txt')
sequence_a = file_a.read()
file_a.close()

file_b = open('fig-6-seq-2.txt')
sequence_b = file_b.read()
file_b.close()

features_a=[
    GraphicFeature(start=495, end=1505, strand=+1, color="tab:olive"),
    GraphicFeature(start=1683, end=1985, strand=+1, color="tab:olive")
]
features_b=[
    GraphicFeature(start=1, end=354, strand=+1, color="tab:olive"),
    GraphicFeature(start=909, end=3824, strand=+1, color="tab:olive"),
    GraphicFeature(start=3887, end=4168, strand=+1, color="tab:olive"),
    GraphicFeature(start=4170, end=4445, strand=+1, color="tab:olive"),
    GraphicFeature(start=4526, end=4819, strand=-1, color="tab:olive")
]
#features[0].label = "ORF1"
#features[0].fontdict["size"] = 9
#features[1].label = "ORF2"
#features[1].fontdict["size"] = 9
#%%
plasmid_a = CircularGraphicRecord(sequence=sequence_a, sequence_length=2378, features=features_a)
plasmid_b = CircularGraphicRecord(sequence=sequence_b, sequence_length=5214, features=features_b)
#fig, _ = plasmid_a.plot(figure_width=4, figure_height=4)
#fig.set_xticks(np.arange(0, 2800, 900))

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(5, 3))#, facecolor="None")
plasmid_a.plot(ax=ax1)
plasmid_b.plot(ax=ax2)
ax1.text(0.5, 0.5, "pCA2.4\n2378 bp", transform=ax1.transAxes, size=9, va='top', ha='center')
ax1.text(1, 0.2, "ORF1", transform=ax1.transAxes, size=9, va='top', ha='right')
ax1.text(0.1, 0.5, "ORF2", transform=ax1.transAxes, size=9, va='top', ha='right')
ax1.text(0, 1, "(a)", transform=ax1.transAxes, size=9, va='top', ha='right')
ax2.text(0.5, 0.5, "pCC5.2\n5214 bp", transform=ax2.transAxes, size=9, va='top', ha='center')
ax2.text(0.55, 0.9, "ORFa", transform=ax2.transAxes, size=9, va='top', ha='center')
ax2.text(0.9, 0.2, "ORFb", transform=ax2.transAxes, size=9, va='top', ha='center')
ax2.text(0.09, 0.52, "ORFc", transform=ax2.transAxes, size=9, va='top', ha='right')
ax2.text(0.12, 0.65, "ORFd", transform=ax2.transAxes, size=9, va='top', ha='right')
ax2.text(0.25, 0.8, "ORFf", transform=ax2.transAxes, size=9, va='top', ha='right')
ax2.text(0, 1, "(b)", transform=ax2.transAxes, size=9, va='top', ha='right')
fig.tight_layout()
#fig.savefig("../figure/thesis-sp.pdf")
fig.savefig("../figure/fig-6.png", dpi=300, bbox_inches='tight', transparent=True)
