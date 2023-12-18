### Script to map identified peptides to BED format
# First read peptide list and proteins to extract from the GFF
# Need a regex to extra protein IDs correctly
# Keep these CDS blocks only
# BED format
# chr_num   chr_start   chr_end name score[0-1000]  strand(+/-) thickStart  thickEnd(same as chr)   0   blockStarts blockSizes  protein peptide unique/not
# Species GenomVer   q-value mod charge  exp_mass   calc_mass   rank    PXD URI

import pandas as pd
from Bio import SeqIO

#Parse fasta and just keep proteins we need
def fasta_parser(fasta_file,prot_pep_dict):
    records = SeqIO.parse(fasta_file, "fasta")
    fasta_id_to_seq = {}

    for record in records:
        if record.id in prot_pep_dict:
            fasta_id_to_seq[record.id] = str(record.seq)
    return fasta_id_to_seq


def gff_parser(gff_file,pep_dict,prot_pep_dict,id_to_protseq,bed_out):

    protein_id_cds_rows = {}
    f = open(gff_file,"r")
    for line in f:
        line = line.rstrip()
        if line.startswith("#") == False:
            cells = line.split("\t")
            if len(cells) > 5:   #not an error
                if cells[2] == "CDS":
                    protein_id = cells[8].split("=")[1] #Parent=Os01t0100100-01
                    #ToDo - GFF3 parser needs work, depending on ID format
                    print("prot:",protein_id)
                    if protein_id in prot_pep_dict:
                        cds_rows = []
                        if protein_id in protein_id_cds_rows:
                            cds_rows = protein_id_cds_rows[protein_id]
                        cds_rows.append(cells)
                        protein_id_cds_rows[protein_id] = cds_rows

    f_out = open(bed_out,"w")
    f_out.write("chr_num\tchr_start\tchr_end\tname\ttscore[0-1000]\tstrand(+/-)\tthickStart\tthickEnd(same as chr)\treserved(0)\t"
                "blockCount\tblockSizes\tchromStarts\tprotein\tpeptide;start;end(removestartendforrealbed)\tunique/not\n")
    pep_counter = 0
    for peptide in pep_dict:
        proteins = pep_dict[peptide]
        for protein in proteins:
            cds_rows = protein_id_cds_rows[protein]
            if protein in id_to_protseq:
                protein_seq = id_to_protseq[protein]
                #f_out.write(make_bed_row(cds_rows,peptide,protein,protein_seq))
                f_out.write(make_bed_row(cds_rows, peptide, protein, protein_seq))
                pep_counter +=1

                if pep_counter ==2000:
                    exit()

            else:
                print("Error, protein seq not found in dictionary2")



### Need to go back to the drawing board;
### Not sure if any of the block positions or block Lengths are correct
### Re-design a simpler algorithm !!



def make_bed_row(cds_list,pep_seq,prot_id,prot_seq):
    first_cds = cds_list[0]
    chr_num = first_cds[0]
    strand = first_cds[6]
    bed_row = ""
    if prot_seq.count(pep_seq) > 1:
        print("Fatal error - haven't yet implemented peptides mapping to multiple positions in one protein")
        exit()
    elif prot_seq.count(pep_seq) == 0:
        print(pep_seq,"not found in",prot_seq)
        exit()

    #Example PEPTYDER Length = 8, mapped to protein positions 2:9  MPEPTYDER....
    #Pepstart =2
    pep_start = prot_seq.find(pep_seq) + 1

    #pep_end = 2 + 8 - 1 = 9; Check
    pep_end = pep_start + len(pep_seq) - 1  # Need to take away one, since a 2 amino acid protein would start at 30 and end at 31, inclusive counting
    print("pep_Start",pep_start,"pep_end",pep_end)
    if strand == "+":

        #Let's assume chromosome start of CDS = 10000#
        #dna_start_offset = (2 * 3) -3 = 3
        #chromStart will be 10000 + 3 = 10003; check
        #dnd_length = 8 * 3 = 24
        #ATGCCCGGGGGGGGGGGGGGGGGGCGA    #Only used proper codons for MPXXXXXXR
        #chromEnd (1 exon case) = chrom_start + dna_length = 100003 + 24 -1 = 100026 #for inclusive counting
        #
        #CCC = P at position
        dna_start_offset = (pep_start * 3) - 3  # Move back 3 to get zero based first base
        #dna_end_offset = (pep_end * 3)  # Aligns already with final base

        for cds in cds_list:
            cds_start = int(cds[3])
            cds_end = int(cds[4])
            cds_len = cds_end - cds_start + 1  # Inclusive count
            #current_cds_relative_end += cds_len

            if dna_start_offset >= cds_start and dna_start_offset <= cds_end:

                #print("chr_start=",chr_start)


    else:
        print("Negative strand not implemented")
    return bed_row



#Will make a single proBED row, given a list of lists of CDS elements, peptide sequence and protein sequence
def make_bed_row2(cds_list,pep_seq,prot_id,prot_seq):
    bed_line = ""
    #First get unchanging metadata
    pep_start = prot_seq.find(pep_seq)  #MPEPTIDERKRKR
    pep_end = pep_start + len(pep_seq)-1    #Need to take away one, since a 2 amino acid protein would start at 30 and end at 31, inclusive counting
    first_cds = cds_list[0]
    chr_num = first_cds[0]
    strand = first_cds[6]
    cds_1_start = 0

    cds_1_start = int(first_cds[3])
    dna_start_offset = (pep_start * 3) - 3  #Move back 3 to get zero based first base
    dna_end_offset = (pep_end * 3)          #Aligns already with final base
    total_blocks_length = dna_end_offset - dna_start_offset

    #pos_tracker = cds_1_start + dna_start_offset

    current_cds_relative_start = 1
    current_cds_relative_end = -1   #To be set
    cds_count = 1
    total_cds_mapped = 1
    start_mapped_cds = 0
    end_mapped_cds = 0
    bed_start = 1
    bed_end = 1

    chr_start = -1
    chr_end = -1
    chrom_starts = []
    blocks_sizes = []

    for cds in cds_list:
        cds_start = int(cds[3])
        cds_end = int(cds[4])
        cds_len = cds_end - cds_start + 1   #Inclusive count
        current_cds_relative_end += cds_len

        #print(cds_start,cds_end,dna_start_offset,dna_end_offset,prot_id,pep_seq)
        
        #
        if dna_start_offset <= current_cds_relative_end and dna_start_offset >= current_cds_relative_start:


            if strand == "+":
                start_within_cds = dna_start_offset - current_cds_relative_start
                chr_start = cds_start + start_within_cds
            else:
                start_within_cds = dna_start_offset - current_cds_relative_end
                chr_start = cds_end - start_within_cds

            block_len_in_cds = cds_len - start_within_cds
            #End reached
            if block_len_in_cds > total_blocks_length:
                block_len_in_cds = total_blocks_length
            chrom_starts.append("0")   #Make string lists since they are only going to be printed
            blocks_sizes.append(str(block_len_in_cds))
        if dna_end_offset <= current_cds_relative_end and dna_end_offset >= current_cds_relative_start:
            if strand == "+":
                end_within_cds = dna_end_offset - current_cds_relative_start
                chr_end = cds_start + end_within_cds
            else:
                end_within_cds = dna_end_offset - current_cds_relative_start
                chr_end = cds_end - end_within_cds  #TODO not really sure this is right, needs testing

            block_start = chr_start - cds_start #This is wrong

            block_len_in_cds = end_within_cds
            if block_len_in_cds < total_blocks_length:  #handle the chunk in this CDS, otherwise all handled already
                blocks_sizes.append(str(end_within_cds))
                chrom_starts.append(str(block_start))
        elif dna_end_offset < current_cds_relative_start and dna_end_offset > current_cds_relative_end: #case for peptides that map more than 2 exons - TODO test this!
            blocks_sizes.append(str(cds_len))
            chrom_starts.append(str(current_cds_relative_start))
        current_cds_relative_start = current_cds_relative_end + 1   #End of CDS block
        cds_count+=1

    bed_line += str(chr_num) + "\t" + str(chr_start) + "\t" + str(chr_end) + "\t" + prot_id+"_"+str(pep_start) + "\t"
    bed_line += "1000" +"\t" + strand + "\t" + str(chr_start) + "\t" + str(chr_end) + "\t"
    bed_line += "0" + "\t" + str(len(chrom_starts)) + "\t" +  str(",".join(blocks_sizes)) + "\t"
    bed_line += str(",".join(chrom_starts)) +  "\t" + prot_id + "\t" + pep_seq + ";" + str(pep_start) + ";" + str(pep_end) +"\t1" + "\n"

    #print(bed_line)
    #exit()
    #print("\t\ttotal blocks=",total_cds_mapped,sep="")



    return bed_line

#Will make a single proBED row, given a list of lists of CDS elements, peptide sequence and protein sequence
def make_bed_row_positive_strand_working(cds_list,pep_seq,prot_id,prot_seq):
    bed_line = ""
    #First get unchanging metadata
    pep_start = prot_seq.find(pep_seq)  #MPEPTIDERKRKR
    pep_end = pep_start + len(pep_seq)-1    #Need to take away one, since a 2 amino acid protein would start at 30 and end at 31, inclusive counting
    first_cds = cds_list[0]
    chr_num = first_cds[0]
    strand = first_cds[6]
    cds_1_start = 0
    if strand == "+":
        cds_1_start = int(first_cds[3])
        dna_start_offset = (pep_start * 3) - 3  #Move back 3 to get zero based first base
        dna_end_offset = (pep_end * 3)          #Aligns already with final base
        total_blocks_length = dna_end_offset - dna_start_offset

        #pos_tracker = cds_1_start + dna_start_offset

        current_cds_relative_start = 1
        current_cds_relative_end = -1   #To be set
        cds_count = 1
        total_cds_mapped = 1
        start_mapped_cds = 0
        end_mapped_cds = 0
        bed_start = 1
        bed_end = 1

        chr_start = -1
        chr_end = -1
        chrom_starts = []
        blocks_lengths = []

        for cds in cds_list:
            cds_start = int(cds[3])
            cds_end = int(cds[4])
            cds_len = cds_end - cds_start + 1   #Inclusive count
            current_cds_relative_end += cds_len

            #print(cds_start,cds_end,dna_start_offset,dna_end_offset,prot_id,pep_seq)
            if dna_start_offset <= current_cds_relative_end and dna_start_offset >= current_cds_relative_start:
                start_within_cds = dna_start_offset - current_cds_relative_start
                chr_start = cds_start + start_within_cds

                block_len_in_cds = cds_len - start_within_cds
                #End reached
                if block_len_in_cds > total_blocks_length:
                    block_len_in_cds = total_blocks_length
                chrom_starts.append("0")   #Make string lists since they are only going to be printed
                blocks_lengths.append(str(block_len_in_cds))
            if dna_end_offset <= current_cds_relative_end and dna_end_offset >= current_cds_relative_start:
                end_within_cds = dna_end_offset - current_cds_relative_start
                chr_end = cds_start + end_within_cds
                block_start = cds_start - chr_start

                block_len_in_cds = end_within_cds
                if block_len_in_cds < total_blocks_length:  #handle the chunk in this CDS, otherwise all handled already
                    blocks_lengths.append(str(end_within_cds))
                    chrom_starts.append(str(block_start))
            elif dna_end_offset < current_cds_relative_start and dna_end_offset > current_cds_relative_end: #case for peptides that map more than 2 exons - TODO test this!
                blocks_lengths.append(str(cds_len))
                chrom_starts.append(str(current_cds_relative_start))
            current_cds_relative_start = current_cds_relative_end + 1   #End of CDS block
            cds_count+=1

        bed_line += str(chr_num) + "\t" + str(chr_start) + "\t" + str(chr_end) + "\t" + prot_id+"_"+str(pep_start) + "\t"
        bed_line += "1000" +"\t" + strand + "\t" + str(chr_start) + "\t" + str(chr_end) + "\t"
        bed_line += "0" + "\t" + str(len(chrom_starts)) + "\t" +  str(",".join(blocks_lengths)) + "\t"
        bed_line += str(",".join(chrom_starts)) +  "\t" + prot_id +  "\t"+pep_seq + "\t" + str(pep_start) + "\t" + str(pep_end) + "\n"

        print(bed_line)
        exit()

    else:
        print("minus strand not yet implemented")
        #TODO
        #TODO have to implement sorting

    return bed_line


### Arguments are the peptide-level tsv file (pep sequence header = "peptide") and the column name containing proteins (semi-colon separated)
def build_pep_prot_maps(pep_tsv_file,protein_fasta_column):

    pep_to_prot = {}
    prot_to_pep = {}
    df = pd.read_csv(pep_tsv_file,sep="\t")

    for i in range(0,len(df)):
        peptide = df.loc[i,"peptide"]
        if pd.isna(df.loc[i,protein_fasta_column]) == False:
            mapped_proteins = df.loc[i,protein_fasta_column].split(";")
            #Make bidirectional dictionaries
            protein_list = []
            if peptide in pep_to_prot:
                protein_list = pep_to_prot[peptide]
            protein_list = protein_list + mapped_proteins
            pep_to_prot[peptide] = protein_list
            for protein in protein_list:
                peptide_list = []
                if protein in prot_to_pep:
                    peptide_list = prot_to_pep[protein]
                peptide_list.append(peptide)
                prot_to_pep[protein] = peptide_list
        else:
            print("No proteins mapped for peptide",peptide)

    return [pep_to_prot,prot_to_pep]



#pep_file = "D:/temp/PXD028712/PXD028712_pep_prophet_thresholded_final_collapsed_pep_level__thresholded_prot_mapped.tsv"
pep_file = "D:/Dropbox/DocStore/Rice/PeptideAtlas/bed_testing/peptide_results.txt"
#ToDo - these results have a final column called RAP-DB, not sure where it came from

#protein_col = "RAPDB"
#gff_file = "D:/Dropbox/DocStore/Rice/PeptideAtlas/bed_testing/Oryza_sativa.IRGSP-1.0.55.chr.gff3"
protein_col = "protein"
gff_file = "D:/Dropbox/DocStore/Rice/PeptideAtlas/bed_testing/mini.gff3"
fasta_location = "D:/Dropbox/DocStore/Rice/PeptideAtlas/bed_testing/fake_protein.fasta"
bed_file = "D:/Dropbox/DocStore/Rice/PeptideAtlas/bed_testing/peptide_results.bed"
#fasta_location = "D:/Dropbox/DocStore/Rice/PanOryza/MAGIC16/fasta/protein/protein/IRGSP-1.0_protein_2022-03-11.fasta"
#bed_file = "D:/temp/PXD028712/PXD028712_pep_prophet_thresholded_final_collapsed_pep_level__thresholded_prot_mapped.bed.tsv"

[pep_2_prot,prot_2_pep] = build_pep_prot_maps(pep_file,protein_col)
id_to_seq_dict = fasta_parser(fasta_location,prot_2_pep)
gff_parser(gff_file,pep_2_prot,prot_2_pep,id_to_seq_dict,bed_file)


