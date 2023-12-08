### Read local CSV format and take best scoring PEPTIDE ID, and collapse - (keeping count of PSMs with Prob > threshold)
### Might want to trial use of pandas? ###

import pandas as pd
import numpy as np
import sys

def read_and_collapse(infile):

    df = pd.read_csv(infile,sep="\t")
    #pivot_df = df.pivot_table(index=["peptide"],values = ["prob"], aggfunc={np.max,np.count_nonzero})
    df_only_cols = df[["mod_peptide","prob"]]
    df_only_cols = df_only_cols.rename(columns={"prob":"count_of_PSMs"})
    group_df = df_only_cols.groupby(by="mod_peptide").count()
    outfile = infile[:-4] + "_count.tsv"
    #group_df.to_csv(outfile,sep="\t")

    #outfile = infile[:-4] + "_collapsed.tsv"
    #pivot_df.to_csv(outfile,sep="\t",index=False)

    df = df.sort_values("prob", ascending=False)    #check sorted, I think it will have been sorted by the Prophets anyway
    df_no_dups = df.drop_duplicates(subset="mod_peptide", keep='first', inplace=False)
    df_no_dups = df_no_dups.set_index("mod_peptide")
    outfile = infile[:-4] + "_no_dups.tsv"
    #df_no_dups.to_csv(outfile,sep="\t")

    final_df = pd.concat([df_no_dups,group_df],axis=1)
    outfile = infile[:-4] + "_final_collapsed.tsv"
    final_df.to_csv(outfile, sep="\t",index_label="modified_peptide")
    print("Finished - output written to",outfile)




    #df['Protein_loc'] = df.groupby('peptide')['peptide'].transform('count')


def testing():
    outfolder = "D:/temp/PXD028712/"
    #infile= outfolder + "PXD028712_pep_prophet_mini3_thresholded.tsv"
    infile= outfolder + "PXD028712_pep_prophet_thresholded.tsv"
    read_and_collapse(infile)

if len(sys.argv)!= 2:
    print("Exit - expected usage\npython ",  sys.argv[0]," inputlocation/PXD00000_interact-ipro.pep_thresholded.tsv \n"
                                                         "Script collapses based on max probability, and does PSM count with"
                                                         "unique modified peptide string being the unit of grouping")
else:
    read_and_collapse(sys.argv[1])







