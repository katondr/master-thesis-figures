import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as opt
import copy
from matplotlib.ticker import ScalarFormatter
from dna_features_viewer import GraphicFeature, GraphicRecord
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.sans-serif"] = "Times New Roman"
plt.rcParams['font.size'] = "9"

text_file = open('fig-5-seq.txt')
sequence = text_file.read()
text_file.close()

features=[
    GraphicFeature(start=216, end=968, strand=-1, color="gold"),
    GraphicFeature(start=974, end=981, strand=-1, color="tab:blue"),
    GraphicFeature(start=1005, end=1015, strand=0, color="tab:olive"),
    GraphicFeature(start=1011, end=1012, strand=-1, color="tab:grey"),
    GraphicFeature(start=1014, end=1019, strand=-1, color="tab:blue"),
    GraphicFeature(start=1041, end=1047, strand=-1, color="tab:blue"),
    GraphicFeature(start=1071, end=1077, strand=0, color="tab:olive"),
    GraphicFeature(start=1115, end=1135, strand=+1, color="tab:olive"),
    GraphicFeature(start=1154, end=1160, strand=+1, color="tab:blue"),
    GraphicFeature(start=1166, end=1167, strand=+1, color="tab:grey"),
    GraphicFeature(start=1174, end=1181, strand=+1, color="tab:blue"),  
    GraphicFeature(start=1187, end=1768, strand=+1, color="gold")  
]
#features[0].label = "$\it{luxR}$"
#features[1].label = "RBS of $\it{luxR}$"
#features[2].label = "putative luxR binding site"
#features[3].label = "transcription start"
#features[4].label = "-10"
#features[5].label = "-35"
#features[6].label = "CAP/cAMP binding site"
#features[7].label = "inversed repeat"
#features[8].label = "-10"
#features[9].label = "transcription start"
#features[10].label = "RBS of $\it{luxI}$"
#features[11].label = "$\it{luxI}$"

#for i in features:
#    i.fontdict["size"] = 9
#features[0].fontdict["size"] = 9

#%%
#file = entrez.fetch_single_file(["AF170104.1"], None, "nuccore", "gb")
#graphic_record = CDSTranslator().translate_record(file)
#promoter = GraphicFeature(start=1116, end=1170, strand=+1, color="tab:grey")
#promoter.label = "$\it{luxI}$ promoter"
#promoter.fontdict["size"] = 9
graphic_record = GraphicRecord(sequence=sequence, sequence_length=8654, features=features)
#m1 = graphic_record.crop((200, 1800))
#cropped_record.features[2].label = "luxCDABE"
#cropped_record.features.append(promoter)
m1 = graphic_record.crop((960, 1200))
m2 = graphic_record.crop((1000, 1080))
#m3 = graphic_record.crop((1110, 1190))
m3 = graphic_record.crop((1110-35, 1190+5))

#m1.features.append(features)

#for item in graphic_record.__dict__["features"]:
#    item.__dict__["label"] = item.__dict__["label"][3] 
#fig, _ = m1.plot(figure_width=6, figure_height=1.25)
#fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(6, 4))#, facecolor="None")
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8, 4))#, facecolor="None")
#graphic_record.plot(ax=axs[0])
m1.plot(ax=axs[0])
m2.plot(ax=axs[1])
m2.plot_sequence(ax=axs[1], fontdict={"size":6})
m3.plot(ax=axs[2])
m3.plot_sequence(ax=axs[2], fontdict={"size":6})
#m2.plot_sequence(ax=axs[2], fontdict={'size':9})

axs[0].fill_between((1000, 1080), +2, -0.5, alpha=0.2, facecolor="tab:blue")
axs[0].fill_between((1110, 1190), +2, -0.5, alpha=0.2, facecolor="tab:blue")
axs[0].set_xticks(np.arange(960, 1201, 60))
axs[1].set_xticks(np.arange(1000, 1081, 20))
axs[2].set_xticks(np.arange(1110, 1191, 20))

axs[0].text(0.48, 0.7, "$\it{luxR}$ promoter region", transform=axs[0].transAxes, size=9, va='top', ha='right')
axs[0].text(0.64, 0.7, "$\it{luxI}$ promoter region", transform=axs[0].transAxes, size=9, va='top', ha='left')
axs[0].text(0, 0.5, "$\it{luxR}$", transform=axs[0].transAxes, size=9, va='top', ha='left')
axs[0].text(1, 0.5, "$\it{luxI}$", transform=axs[0].transAxes, size=9, va='top', ha='right')
axs[1].text(0.92, 0.7, "CAP/cAMP\nbinding site", transform=axs[1].transAxes, size=9, va='top', ha='center')
axs[1].text(0.56, 0.5, "-35", transform=axs[1].transAxes, size=9, va='top', ha='center')
axs[1].text(0.22, 0.75, "-10", transform=axs[1].transAxes, size=9, va='top', ha='center')
axs[1].text(0.15, 0.8, "TS", transform=axs[1].transAxes, size=9, va='top', ha='center')
axs[1].text(0.075, 0.5, "luxR BS", transform=axs[1].transAxes, size=9, va='center', ha='center')
#axs[2].text(1, 0.5, "$\it{luxI}$", transform=axs[2].transAxes, fontsize=9, va='top', ha='right')
#axs[2].text(0.84, 0.5, "putative RBS", transform=axs[2].transAxes, fontsize=9, va='top', ha='center')
#axs[2].text(0.7, 0.5, "TS", transform=axs[2].transAxes, fontsize=9, va='top', ha='center')
#axs[2].text(0.58, 0.5, "-10", transform=axs[2].transAxes, fontsize=9, va='top', ha='center')
#axs[2].text(0.18, 0.5, "luxR binding site", transform=axs[2].transAxes, fontsize=9, va='top', ha='center')
axs[2].text(1, 0.6, "$\it{luxI}$", transform=axs[2].transAxes, size=9, va='top', ha='right')
axs[2].text(0.85, 0.6, "putative RBS", transform=axs[2].transAxes, size=9, va='top', ha='center')
axs[2].text(0.76, 0.6, "TS", transform=axs[2].transAxes, size=9, va='top', ha='center')
axs[2].text(0.67, 0.6, "-10", transform=axs[2].transAxes, size=9, va='top', ha='center')
axs[2].text(0.41, 0.6, "luxR binding site", transform=axs[2].transAxes, size=9, va='top', ha='center')
#axs[1].text(0.55, 0.5, "promoter\nregion", transform=axs[1].transAxes, fontsize=9, va='top', ha='center')
axs[0].text(0, 0.75, "(a)", transform=axs[0].transAxes, size=9, va='top', ha='left')
axs[1].text(0, 0.75, "(b)", transform=axs[1].transAxes, size=9, va='top', ha='left')
#axs[2].text(0, 0.75, "(c)", transform=axs[2].transAxes, fontsize=9, va='top', ha='left')

#print(dir(m2.features))
#print(m1.features)
#print(type(m1.features[0]))
#%%
fig.tight_layout()
#fig.savefig('../figure/promoter-region.pdf', bbox_inches='tight')
fig.savefig('../figure/fig-5.png', dpi=300, bbox_inches='tight', transparent=True)
