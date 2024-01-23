require(rmotifx)
require(ggplot2)
require(ggseqlogo)
library(Biostrings)
library(stringr)
library(tidyverse)
require(gridExtra)
library(grid)
library(miscset)

# # install rmotifx and ggseqlogo
# devtools::install_github('omarwagih/rmotifx')
# devtools::install_github("omarwagih/ggseqlogo")


wd="C:/Users/jtzve/Desktop/Sufo_Tyrosine/04_MotifEnrichment/"
setwd(wd)

# specify the folder and get all FASTA files to make up our sequence library
fasta_path <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/04_MotifEnrichment/data/SwissProt_FASTA_files/" # Update this with your actual file path
fasta_files <- list.files(fasta_path, pattern = "\\.fasta$", full.names = TRUE)

# use Biostrings readAAStringSet on each file to generate the library
sequences <- lapply(fasta_files, function(file) {
  readAAStringSet(file)
})

# flatten to vector to make working with it easier
sequences_vector <- unlist(lapply(sequences, as.character))
#inspect names
head(names(sequences_vector))

# need to extract SwissProt IDs and reassign as names (everything after sp| and bevore the next |)
names(sequences_vector)[1:10] <- sub("sp\\|([^|]+)\\|.*", "\\1", names(sequences_vector))

# build a dataframe that would have peptidoform ID, peptide sequence, protein ID, and relevant protein sequence
BACKGROUND <- tyrosines_containing_background_inSwissProt[,c(1,6)]

BACKGROUND$peptide <- sub("_.*", "", BACKGROUND$peptidoform_id)

# double-check all peptide IDs are in the sequences we have extracted
sum(BACKGROUND$cleaned_protein_IDs %in% names(sequences_vector))
#14718 so yes
names(BACKGROUND)


# extract_surrounding_seq <- function(full_seq, peptide_seq) {
#   match_pos <- gregexpr(peptide_seq, full_seq, fixed = TRUE)[[1]]
#   if (match_pos == -1) {
#     return("not matched") # Return "not matched" if the peptide sequence is not found
#   } else {
#     # Adjust start and end positions to stay within sequence boundaries
#     start <- max(1, match_pos - 7)
#     end <- min(nchar(full_seq), match_pos + nchar(peptide_seq) - 1 + 7)
#     return(substr(full_seq, start, end))
#   }
# }


extract_surrounding_seq <- function(full_seq, peptide_seq) {
  match_positions <- gregexpr(peptide_seq, full_seq, fixed = TRUE)[[1]]
  
  if (match_positions[1] == -1) {
    return("not matched") # Return "not matched" if the peptide sequence is not found
  } else {
    # Handle multiple matches
    surrounding_seqs <- sapply(match_positions, function(pos) {
      start <- max(1, pos - 7)
      end <- min(nchar(full_seq), pos + nchar(peptide_seq) - 1 + 7)
      substr(full_seq, start, end)
    })
    
    # Combine multiple matches into a single string (if needed)
    return(paste(surrounding_seqs, collapse = "; "))
  }
}


# test <- BACKGROUND[1:10,]
# test1 <- test
# test2 <- test
# 
# # Initialize the new column
# test$surrounding_seq <- NA
# 
# for (i in 1:nrow(test)) {
#   peptide_id <- test$cleaned_protein_IDs[i]
#   peptide_seq <- test$peptide[i]
# 
#   if (peptide_id %in% names(sequences_vector)) {
#     full_seq <- sequences_vector[peptide_id]
#     test$surrounding_seq[i] <- extract_surrounding_seq(full_seq, peptide_seq)
#   } else {
#     warning(paste("Peptide ID", peptide_id, "not found in sequences."))
#     test$surrounding_seq[i] <- "not matched"
#   }
# }
# 
# ## this worked and luckily the first sequence was a good testcase for when a 
# # peptide is rtight at the end so it works
# 
# 
BACKGROUND$pep_surrounding_seq <- NA

for (i in 1:nrow(BACKGROUND)) {
  peptide_id <- BACKGROUND$cleaned_protein_IDs[i]
  peptide_seq <- BACKGROUND$peptide[i]

  if (peptide_id %in% names(sequences_vector)) {
    full_seq <- sequences_vector[peptide_id]
    BACKGROUND$pep_surrounding_seq[i] <- extract_surrounding_seq(full_seq, peptide_seq)
  } else {
    warning(paste("Peptide ID", peptide_id, "not found in sequences."))
    BACKGROUND$surrounding_seq[i] <- "not matched"
  }
}
# 
# 

extract_surrounding_seq2 <- function(full_seq, peptide_seq) {
  first_Y_pos <- regexpr("Y", peptide_seq)[[1]] # Position of first Y in peptide
  last_Y_pos <- regexpr("Y", str_rev(peptide_seq), fixed = TRUE)[[1]] # Position of last Y in reversed peptide
  last_Y_pos <- nchar(peptide_seq) - last_Y_pos + 1 # Correct position in original peptide
  
  if (first_Y_pos == -1 || last_Y_pos == -1) {
    return("no Y in peptide") # If no Y found in peptide, return this message
  } else {
    full_seq_start_pos <- gregexpr(peptide_seq, full_seq, fixed = TRUE)[[1]][1] # Start of peptide in full sequence
    # Adjusting start and end positions to include surrounding amino acids, considering sequence boundaries
    start_pos <- max(1, full_seq_start_pos + first_Y_pos - 8)
    end_pos <- min(nchar(full_seq), full_seq_start_pos + last_Y_pos + 6)
    
    return(substr(full_seq, start_pos, end_pos))
  }
}


BACKGROUND$tyr_surrounding_seq <- NA

for (i in 1:nrow(BACKGROUND)) {
  peptide_id <- BACKGROUND$cleaned_protein_IDs[i]
  peptide_seq <- BACKGROUND$peptide[i]
  
  if (peptide_id %in% names(sequences_vector)) {
    full_seq <- sequences_vector[peptide_id]
    BACKGROUND$tyr_surrounding_seq[i] <- extract_surrounding_seq2(full_seq, peptide_seq)
  } else {
    warning(paste("Peptide ID", peptide_id, "not found in sequences."))
    BACKGROUND$surrounding_seq[i] <- "not matched"
  }
}
# 
# 









background_2<- distinct(background)
names(sequences_vector)



######

# Read in sequences
background<-read.csv(file="motif_background.txt", header=FALSE)
bg.seqs = as.character(background$V1)
bg_list<-unique(as.character(background$V3))

file_list<-c("gsb_motif_seqs.txt")
for (f in file_list){
  foreground<-read.csv(file=f,header=FALSE)
  fg.seqs = as.character(foreground$V1)
  #find enriched motifs
  mot = motifx(fg.seqs, bg.seqs, central.res = 'ST', min.seqs = 20, pval.cutoff = 1e-6)

  #motif seqlogo of enriched motif sequences
  pos<-c("-7","-6","-5","-4","-3","-2","-1","p","1","2","3","4","5","6","7")
  break_list<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
  plots<-list()
  plots2<-list()
  counter=1
  counter2=1
  mid=round(length(mot$motif)/2)
  for (val in mot$motif){
    m<-grep(val,fg.seqs,value=TRUE)
    motif_plot=ggplot()+geom_logo(m, method="probability")+scale_x_continuous(labels=pos,breaks=break_list)+ggtitle(val)+theme(legend.position = "none")
    plots[[counter]]<-motif_plot
    counter=counter+1
  }  
  
  p1<-grid.arrange(grobs=plots)
}
  