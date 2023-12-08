library(microbenchmark)
library(parallel)
library(ggplot2)
##### benchmarking getting unique PTMs - probably not really needed 
# as it takes a couple of sec per dataset
benchmark_getting_PTMs <- microbenchmark(
  
  old_code = {
          modifications_withamino <- str_extract_all(df$mod_peptide, "[A-Za-z]\\[[0-9]+\\]") %>%
            unlist() %>% # flatten list to vector
            unique() %>% # get only unique ones
            sort() # sort them
  },
  
  new_code = { 
          modifications_withamino <- unique(sort(unlist(stri_extract_all_regex(df$mod_peptide, "[A-Za-z]\\[[0-9]+\\]"))))
  },
  
  times = 10 # You can adjust the number of times each expression is evaluated
)

print(benchmark_getting_PTMs)
autoplot(benchmark_getting_PTMs)


##### benchgmarking PTM counting function.

# Define the function for PTM counting

df2 <- df
df <- df2[1:5000, ]

benchmark_counting_PTMs <- microbenchmark(
  
  old_code = {
    
    PTM_counts <- apply(df, 1, function(row) count_PTMs(row["mod_peptide"], modifications_withamino))
  },
  
  new_code = { 
    # #this uses stringi and parallel computing
    # count_PTMs <- function(peptide, PTMs_vector) {
    #   sapply(PTMs_vector, function(ptm) stri_count_fixed(peptide, ptm)) 
    # }
    # 
    #this uses stringr and old function
    count_PTMs <- function(peptide, PTMs_vector) {
      # for a vector that contains list of PTMs detected in the data
      sapply(PTMs_vector, function(ptm) str_count(peptide, fixed(ptm))) #count the number of times each PTM appears
    }
    # Setting up parallel processing
    no_cores <- detectCores() - 2  # Leave a couple of cores free for system stability
    cl <- makeCluster(no_cores)
    # Load required libraries on each worker node
    clusterEvalQ(cl, library(stringr)) # or library(stringi)
    # Exporting necessary objects and functions to each worker
    clusterExport(cl, varlist = c("count_PTMs", "modifications_withamino", "df"))
    
    
    # Perform the counting in parallel
    PTM_counts <- parApply(cl, df, 1, function(row) count_PTMs(row["mod_peptide"], modifications_withamino))
    
    # Stop the cluster
    stopCluster(cl)
  },
  
  times = 10 # You can adjust the number of times each expression is evaluated
)
benchmark_counting_PTMs_224276rows <- benchmark_counting_PTMs
# benchmark_counting_PTMs_50000rows <- benchmark_counting_PTMs
# benchmark_counting_PTMs_5000rows <- benchmark_counting_PTMs
# benchmark_counting_PTMs_500rows <- benchmark_counting_PTMs
autoplot(benchmark_counting_PTMs_224276rows)


## benchmarkign number of cores

### this si not a good test - setting up the number of cores is slowish compared to computing small datasets
df2 <- df
df <- df2[1:500, ]

count_PTMs <- function(peptide, PTMs_vector) {
  # for a vector that contains list of PTMs detected in the data
  sapply(PTMs_vector, function(ptm) str_count(peptide, fixed(ptm))) #count the number of times each PTM appears
}
# Setting up parallel processing
no_cores <- 8  # Leave a couple of cores free for system stability
cl <- makeCluster(no_cores)
# Load required libraries on each worker node
clusterEvalQ(cl, library(stringr)) # or library(stringi)
# Exporting necessary objects and functions to each worker
clusterExport(cl, varlist = c("count_PTMs", "modifications_withamino", "df"))


benchmark_counting_PTMs <- microbenchmark(
  
  old_code = {
    
    
    # Perform the counting in parallel
    PTM_counts <- parApply(cl, df, 1, function(row) count_PTMs(row["mod_peptide"], modifications_withamino))
    
   
  },
  
  new_code = { 
  
    # Perform the counting in parallel
    PTM_counts <- parApply(cl, df, 1, function(row) count_PTMs(row["mod_peptide"], modifications_withamino))
  
  },
  
  times = 100 # You can adjust the number of times each expression is evaluated
)

# Stop the cluster
stopCluster(cl)
autoplot(twocores)
autoplot(eightcores)

twocores$expr
