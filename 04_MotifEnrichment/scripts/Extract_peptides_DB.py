import re
from Bio import SeqIO
import collections
import pandas as pd

def makehash():
    return collections.defaultdict(makehash)
dict=makehash()

#Add protein positions for semi-tryptic peptides
def tryptic(protein,peptide,position):
    counter=0
    while counter<=len(peptide):
        peptide_temp=peptide[:counter]
        end=len(peptide)-counter
        peptide_temp_reverse=peptide[end:]
        dict[protein][peptide_temp]=position
        position_temp=(len(peptide)-len(peptide_temp_reverse))+position
        dict[protein][peptide_temp_reverse]=position_temp
        counter+=1

#Dict of potential peptides and sites on proteins, used to find PTM protein positions
def extract_peptide_positions(input):
    counter_test=1
    for seq_record in SeqIO.parse(input, "fasta"):
        if "DECOY" not in seq_record.id and "CONTAM" not in seq_record.id:
            #print(counter_test)
            counter_test+=1
            seq = str(seq_record.seq)
            protein = seq_record.id
            # cut at K or R followed by P
            pep_temp=re.sub(r'(?<=[RK])(?=[^P])','\n', seq, 0, re.DOTALL)
            pep=(pep_temp.split())
            protein_length=0
            pos=[]
            peptides = []
            positions = []
            #start position of peptides
            for p in pep:
                start=protein_length
                protein_length+=len(p)
                pos.append(start)
            #add peptides and start pos to dictionary
            for i in zip(pep,pos):
                peptide=(i[0])
                peptides.append(peptide)
                position=(i[1])
                positions.append(position)
                dict[protein][peptide]=position
                #tryptic(protein,peptide,position)

            # add missed cleavages - concat (max 4 per peptide)
            count=0
            for i in peptides:
                pos2 = count
                count+=1
                if pos2 <= (len(peptides) - 5):
                    missed4=(peptides[pos2] + peptides[pos2+1] + peptides[pos2+2] + peptides[pos2+3]+peptides[pos2+4])
                    missed3=(peptides[pos2] + peptides[pos2+1] + peptides[pos2+2] + peptides[pos2+3])
                    missed2=(peptides[pos2] + peptides[pos2+1] + peptides[pos2+2])
                    missed=(peptides[pos2] + peptides[pos2+1])
                    pos_missed=(positions[pos2])
                    dict[protein][missed4]=pos_missed
                    #tryptic(protein,missed4,pos_missed)
                    dict[protein][missed3]=pos_missed
                    #tryptic(protein,missed3,pos_missed)
                    dict[protein][missed2]=pos_missed
                    #tryptic(protein,missed2,pos_missed)
                    dict[protein][missed]=pos_missed
                    #tryptic(protein,missed,pos_missed)
                if pos2 <= (len(peptides) - 4):
                    missed3=(peptides[pos2] + peptides[pos2+1] + peptides[pos2+2] + peptides[pos2+3])
                    missed2=(peptides[pos2] + peptides[pos2+1] + peptides[pos2+2])
                    missed=(peptides[pos2] + peptides[pos2+1])
                    pos_missed=(positions[pos2])
                    dict[protein][missed3]=pos_missed
                    #tryptic(protein,missed3,pos_missed)
                    dict[protein][missed2]=pos_missed
                    #tryptic(protein,missed2,pos_missed)
                    dict[protein][missed]=pos_missed
                    #tryptic(protein,missed,pos_missed)
                if pos2 <= (len(peptides) - 3):
                    missed2=(peptides[pos2] + peptides[pos2+1] + peptides[pos2+2])
                    missed=(peptides[pos2] + peptides[pos2+1])
                    pos_missed=(positions[pos2])
                    dict[protein][missed2]=pos_missed
                    #tryptic(protein,missed2,pos_missed)
                    dict[protein][missed]=pos_missed
                    #tryptic(protein,missed,pos_missed)
                if pos2 <= (len(peptides)-2):
                    missed=(peptides[pos2] + peptides[pos2+1])
                    pos_missed=(positions[pos2])
                    dict[protein][missed]=pos_missed
                    #tryptic(protein,missed,pos_missed)

            #without N-terminal 'M'
            if seq.startswith('M'):
                pep_temp2 = re.sub(r'(?<=[RK])(?=[^P])', '\n', (seq[1:]), 0, re.DOTALL)
                pep2 = (pep_temp2.split())
                protein_length2 = 0
                pos3 = []
                peptides2 = []
                positions2 = []
                for p in pep2:
                    start2 = protein_length2
                    protein_length2 += len(p)
                    pos3.append(start2)
                for i in zip(pep2, pos3):
                    peptide2 = (i[0])
                    peptides2.append(peptide2)
                    position2 = (i[1])
                    positions2.append(position2)
                    dict[protein][peptide2] = position2+1
                    #tryptic(protein,peptide2,position2)
                count2 = 0
                for i in peptides2:
                    pos4 = count2
                    count2 += 1
                    if pos4 <= (len(peptides2) - 5):
                        missed7 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2] + peptides2[pos4 + 3] +peptides2[pos4 + 4])
                        missed6 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2] + peptides2[pos4 + 3])
                        missed5 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2])
                        missed4 = (peptides2[pos4] + peptides2[pos4 + 1])
                        pos_missed2 = (positions2[pos4])
                        dict[protein][missed7] = pos_missed2+1
                        #tryptic(protein, missed7, pos_missed2)
                        dict[protein][missed6] = pos_missed2+1
                        #tryptic(protein, missed6, pos_missed2)
                        dict[protein][missed5] = pos_missed2+1
                        #(protein, missed5, pos_missed2)
                        dict[protein][missed4] = pos_missed2+1
                        #tryptic(protein, missed4, pos_missed2)
                    if pos4 <= (len(peptides2) - 4):
                        missed6 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2] + peptides2[pos4 + 3])
                        missed5 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2])
                        missed4 = (peptides2[pos4] + peptides2[pos4 + 1])
                        pos_missed2 = (positions2[pos4])
                        dict[protein][missed6] = pos_missed2+1
                        #tryptic(protein, missed6, pos_missed2)
                        dict[protein][missed5] = pos_missed2+1
                        #tryptic(protein, missed5, pos_missed2)
                        dict[protein][missed4] = pos_missed2+1
                        #tryptic(protein, missed4, pos_missed2)
                    if pos4 <= (len(peptides2) - 3):
                        missed5 = (peptides2[pos4] + peptides2[pos4 + 1] + peptides2[pos4 + 2])
                        missed4 = (peptides2[pos4] + peptides2[pos4 + 1])
                        pos_missed2 = (positions2[pos4])
                        dict[protein][missed5] = pos_missed2+1
                        #tryptic(protein, missed5, pos_missed2)
                        dict[protein][missed4] = pos_missed2+1
                        #tryptic(protein, missed4, pos_missed2)
                    if pos4 <= (len(peptides2) - 2):
                        missed4 = (peptides2[pos4] + peptides2[pos4 + 1])
                        pos_missed2 = (positions2[pos4])
                        dict[protein][missed4] = pos_missed2+1
                        #tryptic(protein, missed4, pos_missed2)
    name=input.replace(".fasta", "_Peptides_DB.csv")
    with open(name, "w") as f:
        print("Protein Peptide Pos",file=f)
        #for key, nested in sorted(dict.items()):
            #for subkey, value in sorted(nested.items()):
                #if subkey!="":
                    #print(key, subkey, value, file=f)
    return (dict)