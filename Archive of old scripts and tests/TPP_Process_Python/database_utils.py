
from pyteomics import mass, fasta, parser
from os.path import exists
import pandas as pd

def process_fasta(fasta_file,missed_cleavages,min_seq_len):
    #processes fasta file, digests with trypsin
    #TODO allow for other enzymes

    peptides_to_proteins = {}
    proteins_to_peptides={}

    for description, sequence in fasta.FASTA(fasta_file):

        #ToDo - can easily support any enzyme by changing this parameter
        new_peptides = parser.cleave(sequence, 'trypsin',missed_cleavages)
        unique_peps = {}
        nme_new_peptides = parser.cleave(sequence[1:], 'trypsin',missed_cleavages)  #handle NME
        #print("description",description,"peptides",new_peptides)
        prot_acc = description.split(" ")[0]
        proteins_to_peptides[prot_acc] = new_peptides



        for peptide in new_peptides:
            unique_peps[peptide] = 1
            if len(peptide) >= min_seq_len:
                prot_list = []
                if peptide in peptides_to_proteins:
                    prot_list = peptides_to_proteins[peptide]
                prot_list.append(prot_acc)
                peptides_to_proteins[peptide] = prot_list

        for nme_peptide in nme_new_peptides:
            if nme_peptide not in unique_peps:
                if len(nme_peptide) >= min_seq_len:
                    prot_list = []
                    if nme_peptide in peptides_to_proteins:
                        prot_list = peptides_to_proteins[nme_peptide]
                    prot_list.append(prot_acc)
                    peptides_to_proteins[nme_peptide] = prot_list

    return peptides_to_proteins,proteins_to_peptides


def cache_peps_to_protein_map(outfile,pep_to_prot_dict):

    f_out = open(outfile,"w")
    for pep in pep_to_prot_dict.keys():
        prots = pep_to_prot_dict[pep]
        pep_line = pep + "\t"
        for prot in prots:
            pep_line += prot+";"
        pep_line = pep_line[:-1]    #Remove final ;
        pep_line += "\n"
        f_out.write(pep_line)
    f_out.close()

def build_pep_dict_from_cache(cache_file):

    f = open(cache_file,"r")
    pep_to_prot = {}
    for line in f:
        line = line.rstrip()
        cells = line.split("\t")
        pep = cells[0]
        prots = cells[1].split(";")
        pep_to_prot[pep] = prots
    return pep_to_prot


def setup_pep_dict(fasta):
    cache_location = fasta + "_pep_cached.tsv"
    pep_to_prot_dict = {}
    if exists(cache_location):
        pep_to_prot_dict = build_pep_dict_from_cache(cache_location)
        #cache_peps_to_protein_map(cache_location+"_test", pep_to_prot_dict) #roundtrip test
    else:
        [pep_to_prot_dict,prot_to_pep] = process_fasta(fasta,4,6)
        cache_peps_to_protein_map(cache_location,pep_to_prot_dict)
    return pep_to_prot_dict

#TODO - test round trip; identical dictionary after reading from cache

def add_all_protein_maps_to_df(pep_tsv_file,pep_2_prot,fasta_short_name):

    df = pd.read_csv(pep_tsv_file,sep="\t")
    for i in range(0,len(df)):
        pep = df.loc[i,"peptide"]
        if pep in pep_2_prot:
            df.loc[i,fasta_short_name] = ";".join(pep_2_prot[pep])
        else:
            print("peptide not found in cached database",pep)

    pep_out = pep_tsv_file[:-4] + "_prot_mapped.tsv"
    df.to_csv(pep_out,sep="\t")


fasta_location = "D:/Dropbox/DocStore/Rice/PanOryza/MAGIC16/fasta/protein/protein/IRGSP-1.0_protein_2022-03-11.fasta"
fasta_name = "RAPDB"
#pep_to_prot_dict = setup_pep_dict(fasta_location)
pep_file = "D:/temp/PXD028712/PXD028712_pep_prophet_thresholded_final_collapsed_pep_level__thresholded.tsv"
pep_to_prot = setup_pep_dict(fasta_location)
add_all_protein_maps_to_df(pep_file,pep_to_prot,fasta_name)





