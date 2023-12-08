### Script to map to different rice databases and calculate stats

from pyteomics import mass, fasta, parser

import sys
import os


def process_fasta(fasta_file,missed_cleavages):
    #processes fasta file, digests with trypsin
    #TODO allow for other enzymes

    peptides_to_proteins = {}
    proteins_to_peptides={}
    proteins_to_varieties = {}
    peptides_to_varieties = {}

    for description, sequence in fasta.FASTA(fasta_file):

        #ToDo - can easily support any enzyme by changing this parameter
        new_peptides = parser.cleave(sequence, 'trypsin',missed_cleavages)
        #print("description",description,"peptides",new_peptides)
        prot_acc = description.split(" ")[0]
        variety = map_id_to_variety(prot_acc)

        varieties_for_protein = []
        if prot_acc in proteins_to_varieties:
            varieties_for_protein = proteins_to_varieties[prot_acc]
        varieties_for_protein.append(variety)
        proteins_to_varieties[prot_acc] = varieties_for_protein


        proteins_to_peptides[prot_acc] = new_peptides

        for peptide in new_peptides:
            if len(peptide) > 5:    #Not sure if we need 6mers or lower?
                prot_list = []
                if peptide in peptides_to_proteins:
                    prot_list = peptides_to_proteins[peptide]
                prot_list.append(prot_acc)
                peptides_to_proteins[peptide] = prot_list
    return peptides_to_proteins,proteins_to_peptides,proteins_to_varieties

def process_fasta_with_found_peps(fasta_file,missed_cleavages,peptides_to_datasets):
    #processes fasta file, digests with trypsin
    #TODO allow for other enzymes

    peptides_to_proteins = {}
    proteins_to_peptides={}
    proteins_to_variety = {}
    peptides_to_varieties = {}
    counter = 0
    for description, sequence in fasta.FASTA(fasta_file):

        if counter % 10000 == 0:
            print("Proteins processed:",counter)
        #ToDo - can easily support any enzyme by changing this parameter
        new_peptides = parser.cleave(sequence, 'trypsin',missed_cleavages)
        #print("description",description,"peptides",new_peptides)
        prot_acc = description.split(" ")[0]
        variety = map_id_to_variety(prot_acc)
        proteins_to_variety[prot_acc] = variety

        proteins_to_peptides[prot_acc] = new_peptides

        for peptide in new_peptides:
            if peptide in peptide_to_datasets:    #Not sure if we need 6mers or lower?
                prot_list = []
                if peptide in peptides_to_proteins:
                    prot_list = peptides_to_proteins[peptide]
                prot_list.append(prot_acc)
                peptides_to_proteins[peptide] = prot_list
        counter += 1
    return peptides_to_proteins,proteins_to_peptides,proteins_to_variety


def map_id_to_variety(id):
    variety = ""
    rap_db_set = {"Os01t","Os02t","Os03t","Os04t","Os05t","Os06t","Os07t","Os08t","Os09t","Os10t","Os11t", "Os12t","Os13t"}

    if id[0:6] == "LOC_Os":
        variety = "MSU"
    elif id[0:5] in rap_db_set:
        variety = "RAP"
    elif id[0:4] == "Osaz":
        variety = "Azu"
    elif id[0:4] == "Osat":
        variety = "Nip"
    elif id[0:4] == "Osir":
        variety = "IR64"
    elif id[0:4] == "Osmh":
        variety = "MH63"
    elif id[0:4] == "Oszs":
        variety = "ZS97"
    elif id[0:4] == "Osau":
        variety = "N22"
    elif id[0:8] == "Os117425":
        variety = "ARC"
    elif id[0:8] == "Os127652":
        variety = "NaBo"
    elif id[0:8] == "Os127742":
        variety = "PR106"
    elif id[0:8] == "Os125827":
        variety = "LiXu"
    elif id[0:8] == "Os127564":
        variety = "Lima"
    elif id[0:8] == "Os125619":
        variety = "LaMu"
    elif id[0:8] == "Os127518":
        variety = "KYG"
    elif id[0:8] == "Os128077":
        variety = "KeNa"
    elif id[0:8] == "Os132424":
        variety = "GoSa"
    elif id[0:8] == "Os132278":
        variety = "CMeo"
    elif id[0:5] == "DECOY":
        variety = "decoy"
    elif id[0:5] == "ChrSy":
        variety = "MSU"
    elif id[0:5] == "ChrUn":
        variety = "MSU"
    else:
        print("ID not classified",id)
        variety = "contaminant"
    return variety

if len(sys.argv)!= 4:
    print("Exit - expected usage\npython ",  sys.argv[0],"INPUT_FASTA INPUT_peps_to_prots.txt OUTPUT_TEST\n")
    exit()

fasta_file = sys.argv[1]
in_file = sys.argv[2]
out_file = sys.argv[3]

def run_tests():
    fasta_file = "D:/Data/rice_send_to_peptide_atlas/2022-07-18-decoys-magic16_plus_contaminants.fasta"
    out_file = "D:/temp/peptide_to_variety_map.txt"


#peptides_to_proteins,proteins_to_peptides,proteins_to_varieties = process_fasta(fasta_file,4)

f = open(in_file,"r")
peptide_to_datasets = {}
for line in f:
    line = line.rstrip()
    cells = line.split("\t")
    peptide = cells[0]
    datasets = cells[1]
    peptide_to_datasets[peptide] = datasets

peptides_to_proteins,proteins_to_peptides,proteins_to_variety = process_fasta_with_found_peps(fasta_file,4,peptide_to_datasets)





f_out = open(out_file,"w")

for peptide in peptides_to_proteins:
    proteins = peptides_to_proteins[peptide]
    varieties = set()
    for protein in proteins:
        variety = proteins_to_variety[protein]
        varieties.add(variety)
    datasets = peptide_to_datasets[peptide]
    f_out.write(peptide + "\t" + ";".join(proteins) + "\t" + ";".join(varieties) + "\t" + datasets + "\n")

#Algorithm:
#Write a complex Boolean to assign ID to variety name;

#Osazucena
#Osativa.01G00611 - Nip
#Osir64 - IR64
#Osmh63 - MH63
#Oszs97 - ZS97
#Osaus - n22
#Os117425 - ARC 10497
#Os127652 NATAL BORO
#Os127742 PR106
#Os125827  Liu Xu
#Os127564 Lima
#Os125619 Larha Mugad
#Os127518 Khao Yai
#Os128077 Ketananka
#Os132424 Gobal Sail
#Os132278 Chao Meo