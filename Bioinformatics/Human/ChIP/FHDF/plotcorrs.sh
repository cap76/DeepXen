

#multiBigwigSummary bins -b file1.bw file2.bw -o results.npz


#multiBigwigSummary BED-file -b ../SRR5029389_trimmed_sorted_unique.bw  ../SRR5029390_trimmed_sorted_unique.bw  ../SRR5029391_trimmed_sorted_unique.bw  ../SRR5029392_trimmed_sorted_unique.bw \
#H3K4me3.bw H3K9ac.bw H3K27ac.bw H3K36me3.bw H3K27me3.bw H3K4me2.bw H4K20me1.bw H2A.Z.bw H3K4me1.bw H3K9me3.bw \
#-o results.npz --BED selection.bed

multiBigwigSummary bins -b ../SRR5029389_trimmed_sorted_unique.bw  ../SRR5029390_trimmed_sorted_unique.bw  ../SRR5029391_trimmed_sorted_unique.bw  ../SRR5029392_trimmed_sorted_unique.bw \
H3K4me3.bw H3K9ac.bw H3K27ac.bw H3K36me3.bw H3K27me3.bw H3K4me2.bw H4K20me1.bw H2A.Z.bw H3K4me1.bw H3K9me3.bw \
-o results_all.npz 

plotCorrelation \
-in results_all.npz  \
--corMethod pearson --skipZeros \
--whatToPlot scatterplot \
-o scatterplot_PearsonCorr_bigwigScores_all.pdf   \
--outFileCorMatrix PearsonCorr_bigwigScores_all.tab
