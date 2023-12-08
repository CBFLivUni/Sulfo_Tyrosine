
import sys

if len(sys.argv)!= 3:
    print("Exit - expected usage\npython ",  sys.argv[0],"pep_to_datasets.txt  OUTPUT.txt\n")
    exit()

in_file = sys.argv[1]
out_file = sys.argv[2]
f = open(in_file,"r")

dataset_counter_dict = {}

for line in f:
    line = line.rstrip()
    cells = line.split("\t")
    peptide = cells[0]
    datasets = cells[1].split(";")

    for dataset in datasets:
        dataset_counter = 0
        if dataset in dataset_counter_dict:
            dataset_counter = dataset_counter_dict[dataset]
        dataset_counter+= 1
        dataset_counter_dict[dataset] = dataset_counter

f_out = open(out_file,"w")
f_out.write("DataSet ID\tPeptide count\n")
for dataset in dataset_counter_dict:
    f_out.write(dataset + "\t" + str(dataset_counter_dict[dataset]) + "\n")



