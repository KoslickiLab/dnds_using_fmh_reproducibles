from reportlab.lib import colors
from reportlab.lib.colors import white, lightgrey, darkgreen, darkblue
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib.colors import Color
from Bio.SeqFeature import SeqFeature, SimpleLocation

# Helper function to search for a feature by identifier.
def get_feature(features, gene_id, tags=("locus_tag", "gene", "old_locus_tag", "protein_id")):
    for f in features:
        for tag in tags:
            for x in f.qualifiers.get(tag, []):
                if x.strip() == gene_id.strip():
                    return f
    raise KeyError(f"Gene identifier {gene_id} not found in record.")

# New robust helper function to decide if a gene is linked.
def is_gene_linked(feature, linked_set, tags=("protein_id", "locus_tag", "gene", "old_locus_tag")):
    # Normalize the linked_set
    normalized_linked = {x.strip().lower() for x in linked_set}
    for tag in tags:
        for val in feature.qualifiers.get(tag, []):
            norm_val = val.strip().lower()
            # Check if norm_val is exactly equal, or if one string is a substring of the other.
            for linked in normalized_linked:
                if linked == norm_val or linked in norm_val or norm_val in linked:
                    return True
    return False

def get_gene_name(feature, keys=("protein_id", "locus_tag", "gene")):
    for key in keys:
        if key in feature.qualifiers and feature.qualifiers[key]:
            return feature.qualifiers[key][0].strip()
    return ""


# Read Genome A and Genome B from your flat (GenBank) files.
#record_A = SeqIO.read("JAFRFB010000086.1[1..35797].flat", "genbank")
#record_B = SeqIO.read("JAFQXV010000014.1[1..164250].flat", "genbank")

record_A = SeqIO.read("JAFRFB010000086.1[3155..35412].flat", "genbank")
record_B = SeqIO.read("JAFQXV010000014.1[26121..58711].flat", "genbank")
record_C = SeqIO.read("../other_HGT_analyses/orthologs_in_other_species/CP127222.1[1..862443].flat", "genbank")


# Define the cross-link mappings between Genome A and Genome B.
# Tuple: (score, gene_ID_in_A, gene_ID_in_B)

A_vs_B = [
    (0.0250037, "MBR0371492.1", "MBQ2672915.1"),
    (0.001,     "MBR0371511.1", "MBQ2672934.1"), 
    (0.001,     "MBR0371513.1", "MBQ2672936.1"),
    (0.306843,  "MBR0371514.1", "MBQ2672937.1"),
    (0.001,     "MBR0371515.1", "MBQ2672938.1"), 
    (0.001,     "MBR0371516.1", "MBQ2672939.1"),
    (0.0721069, "MBR0371518.1", "MBQ2672941.1"),
    (0.273767,  "MBR0371519.1", "MBQ2672942.1")
]


B_vs_C = {
    (0, "XYZ", "MBQ2672915.1"),
    (0,     "XYZ", "MBQ2672934.1"), 
    (0,     "XYZ", "MBQ2672936.1"),
    (0,  "XYZ", "MBQ2672937.1"),
    (0,     "XYZ", "MBQ2672938.1"), 
    (0,     "XYZ", "MBQ2672939.1"),
    (0, "XYZ", "MBQ2672941.1"),
    (0,  "XYZ", "MBQ2672942.1")
}

A_vs_C = {
    (0, "MBR0371492.1", "XYZ"),
    (0,     "MBR0371511.1", "XYZ"), 
    (0,     "MBR0371513.1", "XYZ"),
    (0,  "MBR0371514.1", "XYZ"),
    (0,     "MBR0371515.1", "XYZ"), 
    (0,     "MBR0371516.1", "XYZ"),
    (0, "MBR0371518.1", "XYZ"),
    (0,  "MBR0371519.1", "XYZ")
}

"""
## Reverse complement
A_vs_B = [
    (0.0250037, "MBR0371492.1", "MBQ2672915.1"),
    (0.001,     "MBR0371513.1", "MBQ2672936.1"),
    (0.306843,  "MBR0371514.1", "MBQ2672937.1"),
    (0.001,     "MBR0371515.1", "MBQ2672938.1") 
]
"""

"""
## Forward strand
A_vs_B = [
    (0.001,     "MBR0371511.1", "MBQ2672934.1"), 
    (0.001,     "MBR0371516.1", "MBQ2672939.1"),
    (0.0721069, "MBR0371518.1", "MBQ2672941.1"),
    (0.273767,  "MBR0371519.1", "MBQ2672942.1")
]
"""


# Build sets of linked gene IDs.
linked_ids_A = {gene_id_A for (_, gene_id_A, _) in A_vs_B}
linked_ids_B = {gene_id_B for (_, _, gene_id_B) in A_vs_B}


# Create a GenomeDiagram to hold the tracks.
diagram = GenomeDiagram.Diagram("Genome A vs Genome B")

# Create two tracks:
track_A = diagram.new_track(2, name="Methanobrevibacter sp.", greytrack=True, height=0.5, start=0, end=len(record_A))
feature_set_A = track_A.new_set()

track_B = diagram.new_track(1, name="Candidatus Saccharibacteria bacterium", greytrack=True, height=0.5, start=0, end=len(record_B))
feature_set_B = track_B.new_set()

# Add dummy features for cross-links (before adding gene arrow features).
for score, gene_id_A, gene_id_B in A_vs_B:
    # Compute link_color by interpolating from midnightblue to white.
    #link_color = colors.linearlyInterpolatedColor(colors.lightblue, colors.white, 0, 1, score)
    intermediate_color = Color(0.89, 0.89, 0.89)
    link_color = intermediate_color
    border_color = intermediate_color

    feat_A = get_feature(record_A.features, gene_id_A)
    dummy_feat_A = SeqFeature(SimpleLocation(feat_A.location.start, feat_A.location.end, strand=0))
    graphic_feat_A = feature_set_A.add_feature(dummy_feat_A, color=link_color, border=border_color)

    feat_B = get_feature(record_B.features, gene_id_B)
    dummy_feat_B = SeqFeature(SimpleLocation(feat_B.location.start, feat_B.location.end, strand=0))
    graphic_feat_B = feature_set_B.add_feature(dummy_feat_B, color=link_color, border=border_color)

    diagram.cross_track_links.append(CrossLink(graphic_feat_A, graphic_feat_B, link_color, border_color))



label_index = 0
font_size=15
custom_labels_A = {"MBR0371492.1":"urvA","MBR0371511.1":"LCP",
#                    "MBR0371513.1":"YidC/Oxa1 membrane insertase","MBR0371514.1":"rnpA","MBR0371515.1":"rpmH",
                    "MBR0371513.1":f"yidC","MBR0371514.1":"","MBR0371515.1":"rnpA\nrpmH",
                    "MBR0371516.1":"DnaA","MBR0371518.1":"RecF","MBR0371519.1":"RNA methyltransferase"}  # This returns your new custom label.

# Add gene arrow features for Genome A.
for feature in record_A.features:
    if feature.type not in ("gene", "CDS"):
        continue
    # First, get a candidate gene name using our helper function.
    basic_gene_name = get_gene_name(feature)
    # If the feature is linked, override with the custom label from custom_labels_A.
    if is_gene_linked(feature, linked_ids_A):
        gene_name = custom_labels_A.get(basic_gene_name, basic_gene_name)
        gene_color = colors.midnightblue  # your custom color for genes of interest
        if gene_name == "RNA methyltransferase":
            angle=90
        else:
            angle=90
        feature_set_A.add_feature(feature, sigil="BIGARROW", color=gene_color,
                                   label=True, name=gene_name, label_size=font_size,
                                   label_angle=angle, stagger=True, position="start")
    else:
        gene_name = basic_gene_name
        gene_color = colors.lightgrey
        feature_set_A.add_feature(feature, sigil="BOX", color=gene_color,
                                   label=False, name=gene_name, label_size=6,
                                   label_angle=0)

# Add gene arrow features for Genome B.

custom_labels_B = {"MBQ2672915.1":"urvA","MBQ2672934.1":"LCP",
                    "MBQ2672936.1":f"yidC","MBQ2672937.1":"","MBQ2672938.1":"rnpA\nrpmH",
                    "MBQ2672939.1":"DnaA","MBQ2672941.1":"RecF","MBQ2672942.1":"RNA methyltransferase"}  # This returns your new custom label.


# Add gene arrow features for Genome A.
for feature in record_B.features:
    if feature.type not in ("gene", "CDS"):
        continue
    # First, get a candidate gene name using our helper function.
    basic_gene_name = get_gene_name(feature)
    # If the feature is linked, override with the custom label from custom_labels_A.
    if is_gene_linked(feature, linked_ids_B):
        gene_name = custom_labels_B.get(basic_gene_name, basic_gene_name)
        gene_color = colors.midnightblue  # your custom color for genes of interest
        if gene_name == "RNA methyltransferase":
            angle=90
        else:
            angle=90
        feature_set_B.add_feature(feature, sigil="BIGARROW", color=gene_color,
                                   label=True, name=gene_name, label_size=font_size,
                                   label_angle=angle, stagger=True)
    else:
        gene_name = basic_gene_name
        gene_color = colors.lightgrey
        feature_set_B.add_feature(feature, sigil="BOX", color=gene_color,
                                   label=False, name=gene_name, label_size=6,
                                   label_angle=0)

# Draw the diagram in a linear format.
#diagram.draw(format="linear", pagesize="A3", fragments=1, start=0, end=max(len(record_A), len(record_B)))
diagram.draw(format="linear", pagesize="A3", fragments=1, start=-1000, end=max(len(record_A), len(record_B))+1000)

diagram.write("A_vs_B_diagram.pdf", "PDF")
