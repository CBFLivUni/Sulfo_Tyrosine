plot_histograms = function(folder, input_file, wd) {
  # 
  # 
  # folder = "PXD000222/PXD000222/U2OS_PLKinh-2h_TiO2-Diemthyl/comet"
  # input_file = "PXD000222_PXD000222_U2OS_PLKinh-2h_TiO2-Diemthy_interact-ipro.pep_calibrated.tsv"
  # 
  # 
  tsv_file = paste(wd, folder, "/", input_file, sep = "")
  df = read.csv(tsv_file, sep = "\t")
  df = subset(df, df$calibrated_error < 0.1) #Keep only sensible errors, remove odd artefacts from isotope matching
  df =subset(df, df$pp_prob > 0.98)#Only keep high probability matches (even though FDR < 0.01 has already been run)
  hist(df$ppm_error, breaks=150)
  hist(df$calibrated_error, breaks=150)
  hist(df$da_error, breaks=150)
  
  # if the folder has / in it, we need to replace it by _ in the pdf name
  if (grepl("/", folder)){  # check if folder has subfolders
    name <- gsub("/", "_", folder) # if it does, reaplce the slashes in the fole path with _ for the pdf name
    output_pdf = paste(wd, folder, "/", "histograms_for_", name, ".pdf", sep="")
    
  } else {
    
    output_pdf = paste(wd, folder, "/", "histograms_for_", folder, ".pdf", sep="")
    
  }
  
  pdf(output_pdf)
  hist(df$ppm_error, breaks=150)
  hist(df$calibrated_error, breaks=150)
  hist(df$da_error, breaks=150)
  dev.off()
}
