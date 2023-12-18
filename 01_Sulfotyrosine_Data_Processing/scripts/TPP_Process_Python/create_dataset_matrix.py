
import sys

if len(sys.argv)!= 4:
    print("Exit - expected usage\npython ",  sys.argv[0],"INPUT_PEP_TO_PROTEIN_FILE.txt pep_to_data_set.txt OUTPUT.txt\n")
    exit()

in_file = sys.argv[1]
in_data_set_file = sys.argv[1]
out_file = sys.argv[3]
f = open(in_file,"r")

peptide_to_variety_table = {}
for line in f:
    line = line.rstrip()
    cells = line.split("\t")
    peptide = cells[0]
    varieties = cells[2].split(";")
    all_varieties = {'KYG': 0, 'GoSa': 0, 'CMeo': 0, 'IR64': 0, 'MH63': 0, 'LiXu': 0, 'NaBo': 0, 'N22': 0, 'Nip': 0,
                     'PR106': 0, 'Lima': 0, 'MSU': 0, 'LaMu': 0, 'ARC': 0, 'Azu': 0, 'decoy': 0, 'contaminant': 0,
                     'ZS97': 0, 'KeNa': 0, 'RAP':0}
    for variety in varieties:
        all_varieties[variety] = 1
    peptide_to_variety_table[peptide] = all_varieties

    #datasets = cells[3]

all_variety_header  = ['KYG', 'GoSa', 'CMeo', 'IR64', 'MH63', 'LiXu', 'NaBo', 'N22', 'Nip', 'PR106', 'Lima', 'MSU', 'LaMu', 'ARC', 'Azu', 'decoy', 'contaminant', 'ZS97', 'KeNa', 'RAP']
all_variety_header.sort()


f_in_datasets = open(in_data_set_file,"r")

peps_to_dataset_count = {}
for line in f_in_datasets:
    line = line.rstrip()
    cells=line.split("\t")
    peptide = cells[0]
    datasets = cells[1].split(";")
    dataset_count = len(datasets)
    peps_to_dataset_count[peptide] = str(dataset_count)



f_out = open(out_file,"w")
f_out.write("peptide" + "\tdatasetcount\t" + "\t".join(all_variety_header)+"\n")
for peptide in peptide_to_variety_table:
    all_varieties = peptide_to_variety_table[peptide]

    f_out.write(peptide + "\t" + peps_to_dataset_count[peptide] + "\t")
    variety_string = ""
    for variety in all_variety_header:
        pep_count_for_var = all_varieties[variety]
        variety_string += str(pep_count_for_var) + "\t"
    variety_string = variety_string[:-1]
    f_out.write(variety_string + "\n")



