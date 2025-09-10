import pandas as pd

# Load the mapping and shared files
mapping_df = pd.read_excel("idmapping_2025_04_08.xlsx")
shared_df = pd.read_csv("shared_uniprot_output.tsv", sep="\t")

# Rename 'From' to 'UniProt ID' in the mapping file
mapping_df.rename(columns={"From": "UniProt ID"}, inplace=True)

# Merge on the now-matching 'UniProt ID' column
merged_df = pd.merge(
    shared_df,
    mapping_df,
    on="UniProt ID",
    how="left"
)

# Select and rename final columns
final_df = merged_df[[
    "UniProt ID",
    "Protein names",
    "Gene Names",
    "Protein families",
    "Query ID File 1",
    "Query ID File 2"
]].rename(columns={
    "Protein names": "Protein Name",
    "Gene Names": "Gene Name",
    "Protein families": "Protein Family"
})

# Save to output
final_df.to_csv("final_uniprot_mapping.tsv", sep="\t", index=False)
print("✅ Final file saved as final_uniprot_mapping.tsv")
