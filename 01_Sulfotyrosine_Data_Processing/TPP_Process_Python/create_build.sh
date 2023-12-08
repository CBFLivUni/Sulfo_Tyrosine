
###1####
#This reads every PXDXXX_interact-ipro.pep_thresholded.tsv it can find, and just puts into a big dictionary
#Missing step. Ideally, we would want to have some kind of combined scored, and re-rank peptides on that basis; taking into account best PSM score and PSM count
#then re-filter for FDR; without this will be inflated FDR. This can already be assessed later in current workflow (decoys carry forward), but cannot be corrected at present,
#since scores don't carry forward.
#Second argument is output file
python3 compile_all_peptide_ids.py peps_to_datasets.txt


####2#####
#This reads the fasta database, and the found peptide ID list from step 1.
#If there is a match while reading the fasta, it includes a mapping from that peptide to that database;
#This uses a complex conditional statement to figure out the source database, based on rules.

python3 Calculate_stats_per_rice_DB.py /home/arjones/hc-storage/PanOryza/proteomics_fasta/2022-07-18-decoys-magic16_plus_contaminants.fasta peps_to_datasets.txt peptide_to_proteins.txt



###3####
# This step summarises presence / absence per peptide, mapped to each variety
python3 create_dataset_matrix.py peptide_to_proteins.txt peps_to_datasets.txt peptide_to_variety_table.txt

###4####
# generate summary stats / counts per dataset
python3 stats_per_dataset.py peps_to_datasets.txt dataset_summary_stats.txt


