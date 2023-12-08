### Analysis of sulfopeptides
source("plot_hist_fun.R")
wd <- paste(getwd(), "/", sep = "")


plot_histograms(folder = "PXD008952",
                input_file = "PXD008952_interact-ipro-ptm.pep_calibrated.tsv",
                wd = wd) #From PeptideAtlas

plot_histograms("CPTAC_S013", "CPTAC_S013_interact-ipro-ptm.pep_calibrated.tsv", wd)
plot_histograms("PXD018273", "PXD018273_interact-prob.pep_calibrated.tsv", wd)   #Zirconium enriched, have only searched small number of Thermo data, there is a large number of Sciex Triple TOF data
plot_histograms("PXD037549", "PXD037549_interact-prob.pep_calibrated.tsv", wd) #Novel phosphospeptide enrichment method
plot_histograms("PXD036069", "PXD036069_interact-prob.pep_calibrated.tsv", wd) #PXD036069 secretome and phospho enriched (mouse) - currently being searched on cluster 
plot_histograms("PXD040897", "PXD040897_interact-prob.pep_calibrated.tsv", wd)


# additional plots from prev datasets
plot_histograms("PXD000222/PXD000222/U2OS_PLKinh-2h_TiO2-Diemthyl/comet", 
                "PXD000222_PXD000222_U2OS_PLKinh-2h_TiO2-Diemthy_interact-ipro.pep_calibrated.tsv", 
                wd)


plot_histograms("PXD000222/PXD000222/U2OS_PLKinh-2h-rep_TiO2-Diemthyl/comet", 
                "PXD000222_PXD000222_U2OS_PLKinh-2h-rep_TiO2-Diemthy_interact-ipro.pep_calibrated.tsv", 
                wd)

plot_histograms("PXD000222/PXD000222/U2OS_PLKinh-4h_TiO2-Diemthyl/comet", 
                "PXD000222_PXD000222_U2OS_PLKinh-4h_TiO2-Diemthy_interact-ipro.pep_calibrated.tsv", 
                wd)

plot_histograms("PXD000222/PXD000222/U2OS_PLKinh-6h_TiO2-Diemthyl/comet", 
                "PXD000222_PXD000222_U2OS_PLKinh-6h_TiO2-Diemthy_interact-ipro.pep_calibrated.tsv", 
                wd)

plot_histograms("PXD000222/PXD000222/U2OS_PLKinh-30m_TiO2-Diemthyl/comet", 
                "PXD000222_PXD000222_U2OS_PLKinh-30m_TiO2-Diemthy_interact-ipro.pep_calibrated.tsv", 
                wd)





