

# Use prodigal to predict protein-coding sequences from scaffolds

prodigal -i scaffold_GCA_017406775.1.fna -o scaffold.GCA_017406775.genes -a scaffold_GCA_017406775.proteins.faa -d scaffold_GCA_017406775.cds.fna 

prodigal -i scaffold_GCA_017434555.1.fna -o scaffold.GCA_017434555.genes -a scaffold_GCA_017434555.proteins.faa -d scaffold_GCA_017434555.cds.fna 

# Use proteinortho to find what protein sequences are shared between scaffolds

proteinortho -project=test -p=blastp+ *faa

# Use diamond

## Build database for swissprot

wget http://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

gunzip uniprot_sprot.fasta.gz

diamond makedb --in uniprot_sprot.fasta -d swissprot


## Running diamond

./diamond blastp -d swissprot -q queries.fasta -o scaffold_GCA_017434555.proteins.matches.tsv

./diamond blastp -d swissprot -q queries.fasta -o scaffold_GCA_017406775.proteins.matches.tsv

# Obtain UNIprot proteins shared between genomes A and B

## Install packages

conda create -n uniprot_env python=3.10
conda activate uniprot_env

conda install bioconda::bioservices
conda install anaconda::pandas 

## Run in-house script to identify shared uniprot protein-coding sequences between Genome A and B

python find_shared_uniprot.py

# Identify uniprot ID protein fmailies

IGNORE THIS: cat shared_uniprot_output.tsv | cut -f 2,4,5 | sort | uniq -c | sort | sed 's/^\s*//' | sed 's/\s/\t/' > sorted_and_counted.tsv

DO THIS INSTEAD: Copy the shared uniprot IDs reported by find_shared_uniprot.py and search under ID mapping in uniprot website, download table, make sure protein family column is included

# Match the uniprot IDs to the scaffolds

python uniprot_info_of_shared_genes.py

cat final_uniprot_mapping.tsv | cut -d $'\t' -f 4,5,6 | sort -t $'\t' -k2,2 -k3,3 -k1,1 | cut -d $'\t' -f 2,3 | sort | uniq -c | sort -n | sed 's/^\s*//' | sed 's/\s/\t/' > total_orthologs_between_genomes.tab

# align, create AXT files, and run kaks_calculator

## Using clustalw align pairs
python align_orthologs.py

## convert alignments to AXT files and run it through kaks_calculator
python convert_alignments_to_axt.py

# produce figure

## figure is a subplot of 3 for Ka, ks, and ka/ks
 
kaks_subplots.pdf