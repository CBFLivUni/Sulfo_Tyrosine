import sys
import statistics

def summary_peptidoform_stats(results_file):
    f = open(results_file,"r")
    rawfile_to_Da_shifts = {}

    counter =0

    peptidoform_col = 0
    error_col = 0
    peptidoform_to_errors = {}
    peptidoform_to_protein = {}
    protein_col = 0
    for line in f:
        line = line.rstrip("\n")
        cells = line.split("\t")
        if counter != 0:
            peptidoform = cells[peptidoform_col]
            protein = cells[protein_col]
            errors = []
            if peptidoform in peptidoform_to_errors:
                errors = peptidoform_to_errors[peptidoform]
            errors.append(float(cells[error_col]))
            peptidoform_to_errors[peptidoform] = errors
            peptidoform_to_protein[peptidoform] = protein

        else:
            for i in range(0,len(cells)):
                if cells[i] == "protein":
                    protein_col = i
                if cells[i] == "calibrated_error":
                    error_col = i
                if cells[i] == "mod_peptide":
                    peptidoform_col = i

        counter+=1

    f_out = open(results_file[:-4]+"_collapse.tsv","w")
    f_out.write("peptidoform\tprotein\tmean_error\tmedian_error\tstdev_error\tcount_peptidoform\tall_errors\n")
    for peptidoform in peptidoform_to_errors:
        errors = peptidoform_to_errors[peptidoform]
        if len(errors) > 2:
            mean_error = str(statistics.mean(errors))
            median_error = str(statistics.median(errors))
            std_dev = str(statistics.stdev(errors))

            f_out.write(peptidoform + "\t")
            f_out.write(peptidoform_to_protein[peptidoform] + "\t")
            f_out.write(mean_error  + "\t")
            f_out.write(median_error + "\t")
            f_out.write(std_dev + "\t")
            f_out.write(str(len(errors)) + "\t")
            f_out.write(";".join([str(ele) for ele in errors]) + "\n")







if len(sys.argv)!= 2:
    print("Exit - expected usage\npython ",  sys.argv[0]," results.tsv\n")
else:
    summary_peptidoform_stats(sys.argv[1])