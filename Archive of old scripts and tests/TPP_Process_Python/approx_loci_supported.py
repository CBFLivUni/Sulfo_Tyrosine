
f = open("testing_pep_index_found_unique_var.txt","r")

msu_proteins = set()
nip_proteins = set()

for line in f:
    line = line.rstrip()
    cells = line.split("\t")
    proteins = cells[1].split(";")

    for protein in proteins:
        if "LOC_Os" in protein:
            gene = protein[:-2] #remove suffix
            msu_proteins.add(gene)
        if "Osativa" in protein:
            gene = protein[:-3]
            nip_proteins.add(gene)

print("msu supported loci:",len(msu_proteins))
print("Nip supported loci:",len(nip_proteins))
