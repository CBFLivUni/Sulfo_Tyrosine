import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import os

all_seq_fg=[]
all_seq_prot_fg=[]
all_motif_bg=[]
all_motif_prot_bg=[]
all_phospho_bg=[]

phospho_sites="D:/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/Rice_phosphosite_matrix_binomial_peptidoform_0.05_w_protein-pos_scores_FLR_SNP_reps_expand.csv"
database="D:/Dropbox/PTMExchange/Rice/database/protein_set_comparisons/protein_set_comparisons/Osativa_super_annotation_union_noIC4R_v2_cRAP.fasta"


df = pd.read_csv(phospho_sites)
#list for all proteins, remove duplicated protein phosphosites
df=df.dropna(subset = ['PTM_pos'], inplace = False)
df=df.loc[df['UP']!="-"]
df['PTM_pos']=df['PTM_pos'].astype(int)
PTM_pos_list=df['PTM_pos'].tolist()
protein_list=df['Protein'].to_list()
category_list=df['PTM_FLR_category'].to_list()

seq_dict=SeqIO.to_dict(SeqIO.parse(database,"fasta",alphabet=IUPAC.extended_protein))
seq_list=[]
seq_protein=[]
# bronze_seq_list=[]
# bronze_seq_protein=[]
# silver_seq_list=[]
# silver_seq_protein=[]
# gold_seq_list=[]
# gold_seq_protein=[]
bsg_seq_list=[]
bsg_seq_protein=[]
sg_seq_list=[]
sg_seq_protein=[]
g_seq_list=[]
g_seq_protein=[]

#All STY phosphosites +/- 7 - 15mer peptides (motif foreground)
for a,b,c in zip(PTM_pos_list,protein_list,category_list):
    if b in seq_dict:
        record=seq_dict[b]
    elif "sp|"+b in seq_dict:
        record=seq_dict["sp|"+b]
    elif "tr|"+b in seq_dict:
        record=seq_dict["tr|"+b]
    else:
        #print(b) #some proteins not present in search DB as using mapped DB from entire MSU DB
        continue
    seq_temp=str(record.seq)


    #remove pA hits
    if seq_temp[a-1]=="A":
        continue

    if a+9>len(seq_temp):
        while a+8>len(seq_temp):
            seq_temp+="_"
    if a<9:
        seq_temp=("_"*(8-a))+seq_temp
        a+=(8-a)

    seq_list.append(seq_temp[a-8:a+7])
    all_seq_fg.append(seq_temp[a-8:a+7])
    seq_protein.append(b)

    #print(seq_temp[a-1])
    #print(seq_temp[a-8:a+7])
    if seq_temp[a-1]!="S" and seq_temp[a-1]!="T" and seq_temp[a-1]!="Y":
        print(a,b,c)
        print(seq_temp)
        print(seq_temp[a-8:a+7])

    if c=="Bronze":
        #bronze_seq_list.append(seq_temp[a-8:a+7])
        #bronze_seq_protein.append(b)
        bsg_seq_list.append(seq_temp[a-8:a+7])
        bsg_seq_protein.append(b)
    if c=="Silver":
        #silver_seq_list.append(seq_temp[a-8:a+7])
        #silver_seq_protein.append(b)
        bsg_seq_list.append(seq_temp[a - 8:a + 7])
        bsg_seq_protein.append(b)
        sg_seq_list.append(seq_temp[a - 8:a + 7])
        sg_seq_protein.append(b)
    if c=="Gold":
        #gold_seq_list.append(seq_temp[a-8:a+7])
        #gold_seq_protein.append(b)
        bsg_seq_list.append(seq_temp[a - 8:a + 7])
        bsg_seq_protein.append(b)
        sg_seq_list.append(seq_temp[a - 8:a + 7])
        sg_seq_protein.append(b)
        g_seq_list.append(seq_temp[a - 8:a + 7])
        g_seq_protein.append(b)

#save as list of sequences - look up Uniprot (DAVID foreground)
UP_dict=dict(zip(df.Protein, df.UP))
#all_seq_list=[seq_list,bronze_seq_list,silver_seq_list,gold_seq_list]
#all_proteins=[seq_protein,bronze_seq_protein,silver_seq_protein,gold_seq_protein]
#names=["motif_seqs.txt","bronze_motif_seqs.txt","silver_motif_seqs.txt","gold_motif_seqs.txt"]
names=["gsb_motif_seqs.txt","gs_motif_seqs.txt","gold_motif_seqs.txt"]
all_seq_list=[bsg_seq_list,sg_seq_list,g_seq_list]
all_proteins=[bsg_seq_protein,sg_seq_protein,g_seq_protein]
for i,j,k in zip(all_seq_list,all_proteins,names):
    df2 = pd.DataFrame(list(zip(i,j)), columns=['Sequences','Proteins'])
    #print(df2)
    updated_protein_list=[]
    for a in range(len(df2)):
        protein=df2.loc[a,'Proteins']
        protein_2=UP_dict[protein]
        #if protein isn't uniprot accession ("|")
        if ("|") in protein_2:
            protein_2=UP_dict[protein].split("|")[1]
        updated_protein_list.append(protein_2)
        if k=="motif_seqs.txt" or k=="gsb_motif_seqs.txt":
            all_seq_prot_fg.append(protein_2)

    df2['Protein']=updated_protein_list
    df2=df2.drop_duplicates()
    df2.to_csv("D:/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/All_datasets/"+k,index=False, header=False)

#print("HERE")
#All STY sites +/- 7 - 15mer (motif background)
background_list=[]
background_protein=[]

for record in SeqIO.parse(database,"fasta",alphabet=IUPAC.extended_protein):
    #print(record)
    seq_temp=record.seq
    for a in range(len(seq_temp)):
        #print(a)
        if seq_temp[a]=="S" or seq_temp[a]=="T" or seq_temp[a]=="S":
            if a+8>len(seq_temp):
                while a+8>len(seq_temp):
                    seq_temp+="_"
            if a<7:
                seq_temp=("_"*(7-a))+seq_temp
                a+=(7-a)
            #print(seq_temp[a-7:a+8])
            background_list.append(seq_temp[a-7:a+8])
            all_motif_bg.append(seq_temp[a-7:a+8])
            background_protein.append(record.id)
            all_motif_prot_bg.append(a)

df2 = pd.DataFrame(list(zip(background_list,background_protein)), columns=['Sequences','Proteins'])
updated_protein_list=[]
for a in range(len(df2)):
    protein=df2.loc[a,'Proteins']
    #print(protein)
    #if protein isn't uniprot accession ("|")
    try:
        protein_2=UP_dict[protein].split("|")[1]
    except:
        protein_2="-"
    #print(protein,protein_2)
    updated_protein_list.append(protein_2)

df4 = pd.DataFrame(list(zip(background_list,background_protein)), columns=['Sequences','Proteins'])
df4['Protein2']=updated_protein_list
#print(len(df4))
df4=df4.loc[df4['Protein2']!="-"]
#print(len(df4))
df4=df4.drop_duplicates()
df4.to_csv("D:/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/All_datasets/motif_background.txt",index=False, header=False)

#All phospho-proteins background
updated_protein_list=[]
for protein in (protein_list):
    #if protein isn't uniprot accession ("|")
    protein_2=UP_dict[protein].split("|")[1]
    updated_protein_list.append(protein_2)
    all_phospho_bg.append(protein_2)
df5=pd.DataFrame(list(updated_protein_list),columns=['Proteins'])
df5=df5.drop_duplicates()
df5.to_csv("D:/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/All_datasets/All_phospoproteins_background.txt",index=False, header=False)

#df6=pd.DataFrame(list(zip(all_seq_fg,all_seq_prot_fg)), columns=['Sequences','Proteins'])
#df6=df6.drop_duplicates()
#df6.to_csv("C:/users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/05.05.23/All_datasets/motif_seqs.txt",index=False, header=False)
#df7 = pd.DataFrame(list(zip(all_motif_bg,all_motif_prot_bg)), columns=['Sequences','Proteins'])
#df7.to_csv("C:/users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated/All_datasets/motif_background.txt",index=False, header=False)
#df8=pd.DataFrame(list(all_phospho_bg),columns=['Proteins'])
#df8=df8.drop_duplicates()
#df8.to_csv("C:/users/krams/Dropbox/PTMExchange/Rice/New_build_ID/FLR_updated//05.05.23/All_datasets/All_phospoproteins_background.txt",index=False, header=False)