### Analysis of sulfopeptides
source("plot_hist_fun.R")
wd <- paste(getwd(), "/", sep = "")


plot_histograms(folder = "PXD008952",
                input_file = "PXD008952_interact-ipro-ptm.pep_calibrated.tsv",
                wd = wd) #From PeptideAtlas

plot_histograms("CPTAC_S013","CPTAC_S013_interact-ipro-ptm.pep_calibrated.tsv", wd)
plot_histograms("PXD018273","PXD018273_interact-prob.pep_calibrated.tsv")   #Zirconium enriched, have only searched small number of Thermo data, there is a large number of Sciex Triple TOF data
plot_histograms("PXD037549","PXD037549_interact-prob.pep_calibrated.tsv") #Novel phosphospeptide enrichment method
plot_histograms("PXD036069","PXD036069_interact-prob.pep_calibrated.tsv") #PXD036069 secretome and phospho enriched (mouse) - currently being searched on cluster 

#plot_histograms("PXD012188","PXD012188_interact-prob.pep_calibrated.tsv") #PXD012188 is Hardman et al EMBO - Claire's non-canonical phospho data - didn't use metal ions to enrich; 29/11 - post-processing on server, has some problems



#PXD040897 Large-scale modern phosphoproteome (Covid infection of host cells)


