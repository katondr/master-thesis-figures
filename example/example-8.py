import matplotlib.pyplot as plt

#font setting
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.sans-serif"] = "Times New Roman"
plt.rcParams['font.size'] = "9"

import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.io.fasta as fasta
import biotite.database.entrez as entrez
import biotite.sequence.graphics as graphics
import biotite.database.entrez as entrez

# Download and parse protein sequences of avidin and streptavidin
fasta_file = fasta.FastaFile.read(entrez.fetch_single_file(
    ["BAA10098", "AAG05773"], None, "protein", "fasta"
))

for name, sequence in fasta_file.items():
    if "BAA10098" in name:
        syn_seq = seq.ProteinSequence(sequence)
    elif "AAG05773" in name:
        pseae_seq = seq.ProteinSequence(sequence)
        
# Get BLOSUM62 matrix
matrix = align.SubstitutionMatrix.std_protein_matrix()
# Perform pairwise sequence alignment with affine gap penalty
# Terminal gaps are not penalized
alignments = align.align_optimal(syn_seq, pseae_seq, matrix,
                                 gap_penalty=(-10, -1), terminal_penalty=False)

# Draw first and only alignment
# The color intensity indicates the similiarity
fig = plt.figure(figsize=(6, 8))
ax = fig.add_subplot(111)
graphics.plot_alignment_similarity_based(
    ax, alignments[0], matrix=matrix, labels=["slr0378", "pvdQ"],
    show_numbers=True, show_line_position=True
)
fig.tight_layout()

#plt.show()
#fig.savefig("../figure/aln-2.pdf", bbox_inches='tight')
fig.savefig("figure/fig-8.png", dpi=300, bbox_inches='tight', transparent=True)
