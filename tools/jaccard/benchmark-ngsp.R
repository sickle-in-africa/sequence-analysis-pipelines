#
#  BENCHMARK R - plotting script
#
#	* plot the output of the benchmark test-suite
#	* as violin plots using ggplot2
#
#	Jack Morrice
#
######################################################
#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
jac_dir = args[2]
plot_prefix = args[3]
multiplot_src = paste(jac_dir,"/multiplot.R", sep="")

# SADaCC colors
# c(red, dark grey, warm grey, brown, warm grey light, light brown)
colors=c("#c62e3a", "#6d6e71", "#c3b9b8", "#cfb9a6", "#d8d7d3", "#e7ceb9")

source(multiplot_src)

df <- read.csv(args[1], header=TRUE, sep="\t")

gg1 <- ggplot(df, aes(aligner, jaccard, color=caller)) +
	geom_violin() + 
	geom_jitter(size=1, shape=16, position=position_jitter(0.1)) +
	labs(y="Jaccard index")

plot_output_fname = paste(plot_prefix, ".jacc.pdf", sep="")
pdf(plot_output_fname)
print(gg1)

gg2 <- ggplot(df, aes(aligner, runtime, color=caller)) +
	geom_violin() + 
	geom_jitter(size=1, shape=16, position=position_jitter(0.1)) +
	labs(y="runtime (s)")

plot_output_fname = paste(plot_prefix, ".time.pdf", sep="")
pdf(plot_output_fname)
print(gg2)

dev.off()