library('diveRsity')

setwd("path/to/wd")

divMigrate(infile = "neutral_pp_snps.gen", outfile = "divMig_0_neutral", boots = 1000, stat = "all",
           filter_threshold = 0, plot_network = TRUE, plot_col = "darkblue", para = TRUE)

divMigrate(infile = "neutral_pp_snps.gen", outfile = "divMig_25_neutral", boots = 1000, stat = "all",
           filter_threshold = .25, plot_network = TRUE, plot_col = "darkblue", para = TRUE)

divMigrate(infile = "neutral_pp_snps.gen", outfile = "divMig_40_neutral", boots = 1000, stat = "all",
           filter_threshold = .40, plot_network = TRUE, plot_col = "darkblue", para = TRUE)

divMigrate(infile = "neutral_pp_snps.gen", outfile = "divMig_35_neutral", boots = 1000, stat = "all",
           filter_threshold = .35, plot_network = TRUE, plot_col = "darkblue", para = TRUE)




