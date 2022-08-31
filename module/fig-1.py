import matplotlib.pyplot as plt
import numpy as np

#plt.rcParams["font.family"] = "serif"
#plt.rcParams["font.sans-serif"] = "Times New Roman"
#plt.rcParams['font.size'] = "9"

from dna_features_viewer import BiopythonTranslator
import biotite.database.entrez as entrez
from dna_features_viewer import compute_features_levels, GraphicRecord, GraphicFeature

class CDSTranslator(BiopythonTranslator):
    """
    """
    
    def compute_feature_fontdict(self, feature):
        """"""
        font = {'family': 'sans-serif', 'weight': 'normal', 'size': 9, 'style': 'italic'}
        return font
    
    def compute_feature_color(self, feature):
        if feature.type == "CDS":
            return "gold"

    def compute_feature_label(self, feature):
        if feature.type == 'CDS' and BiopythonTranslator.compute_feature_label(self, feature) != "unknown":
            return BiopythonTranslator.compute_feature_label(self, feature)
        else:
            return None
    
    def compute_filtered_features(self, features):
        """"""
        return [
            feature for feature in features
            if (feature.type == "CDS") and BiopythonTranslator.compute_feature_label(self, feature) != "unknown"
            #and ("lux" in feature.ref)
        ]

features=[
    GraphicFeature(start=974, end=981, strand=-1, color="tab:grey"),
    GraphicFeature(start=1005, end=1015, strand=0, color="tab:grey"),
    GraphicFeature(start=1011, end=1012, strand=-1, color="tab:grey"),
    GraphicFeature(start=1014, end=1019, strand=-1, color="tab:grey"),
    GraphicFeature(start=1041, end=1047, strand=-1, color="tab:grey"),
    GraphicFeature(start=1071, end=1077, strand=0, color="tab:grey"),
    GraphicFeature(start=1116, end=1135, strand=+1, color="tab:grey"),
    GraphicFeature(start=1154, end=1160, strand=+1, color="tab:grey"),
    GraphicFeature(start=1166, end=1167, strand=+1, color="tab:grey"),
    GraphicFeature(start=1174, end=1181, strand=+1, color="tab:grey"),  
]
features[0].label = "RBS of luxR"
features[1].label = "putative luxR binding site"
features[2].label = "transcription start"
features[3].label = "-10"
features[4].label = "-35"
features[5].label = "CAP/cAMP binding site"
features[6].label = "inversed repeat"
features[7].label = "-10"
features[8].label = "RBS of luxI"

features[0].fontdict["fontsize"] = 9

file = entrez.fetch_single_file(["AF170104.1"], None, "nuccore", "gb")
graphic_record = CDSTranslator().translate_record(file)
#promoter = GraphicFeature(start=1116, end=1170, strand=+1, color="tab:grey")
#promoter.label = "$\it{luxI}$ promoter"
#promoter.fontdict["fontsize"] = 9
    
m1 = graphic_record.crop((200, 1800))
#cropped_record.features[2].label = "luxCDABE"
#cropped_record.features.append(promoter)
m2 = graphic_record.crop((9010, 1250))
#m2.features

for item in graphic_record.__dict__["features"]:
    item.__dict__["label"] = item.__dict__["label"][3] 
    
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6, 3), facecolor="None")
graphic_record.plot(ax=axs[0])
m1.plot(ax=axs[1])
#m2.plot(ax=axs[2])
#m2.plot_sequence(ax=axs[2], fontdict={'size':9})

axs[0].fill_between((215, 1768), +2, -0.5, alpha=0.2, facecolor="tab:blue")
axs[1].fill_between((968, 1187), +2, -0.5, alpha=0.2, facecolor="tab:blue")
axs[0].set_xticks(np.arange(0, 9000, 2000))
#axs[1].set_xticks(np.arange(0, 2100, 500))
axs[0].text(1, 0.75, "sequence length: "+str(graphic_record.sequence_length), transform=axs[0].transAxes, fontsize=9, va='top', ha='right')
axs[0].text(0.12, 0.5, "regulatory\ngenes", transform=axs[0].transAxes, fontsize=9, va='top', ha='center')
axs[1].text(0.55, 0.5, "promoter\nregion", transform=axs[1].transAxes, fontsize=9, va='top', ha='center')
axs[0].text(-0.025, 0.5, "(a)", transform=axs[0].transAxes, fontsize=9, va='top', ha='right')
axs[1].text(-0.025, 0.5, "(b)", transform=axs[1].transAxes, fontsize=9, va='top', ha='right')

fig.tight_layout()
#fig.savefig('../figure/lux-operon-map.pdf', bbox_inches='tight')
fig.savefig('../figure/lux-operon-map.png', dpi=300, bbox_inches='tight', transparent=True)
