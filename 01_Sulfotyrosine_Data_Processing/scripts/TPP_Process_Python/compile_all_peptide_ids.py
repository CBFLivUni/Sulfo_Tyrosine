import os

import sys

if len(sys.argv)!= 2:
    print("Exit - expected usage\npython ",  sys.argv[0],"peps_to_datasets.txt\n")
    exit()

outfile = sys.argv[1]

folders = []
for file in os.listdir("."):
    if os.path.isdir(file):
        folders.append(file)

peptide_to_datasets = {}
for f in folders:
    #PXD031352_interact-ipro.pep.tsv
    results_file = f + "/" + f + "_interact-ipro.pep_thresholded.tsv"
    dataset_id = f

    if os.path.isfile(results_file):
        print("processing",results_file)
        f = open(results_file,"r")
        counter = 0
        for line in f:
            line = line[:-1]
            if counter != 0:
                cells = line.split("\t")
                peptide = cells[6]
                #print(peptide)
                datasets_for_pep = set()
                if peptide in peptide_to_datasets:
                    datasets_for_pep = peptide_to_datasets[peptide]
                datasets_for_pep.add(dataset_id)
                peptide_to_datasets[peptide]=datasets_for_pep

            counter +=1
        print("\t" + results_file + " = ",counter," PSMs")
        print("\ttotal peptides in dict:", len(peptide_to_datasets))

f_out = open(outfile,"w")
print("num peptides",len(peptide_to_datasets))
for peptide in peptide_to_datasets:
    datasets_for_pep = peptide_to_datasets[peptide]
    f_out.write(peptide + "\t" + ";".join(datasets_for_pep) + "\n")
f_out.close()
