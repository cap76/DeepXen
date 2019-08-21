computeMatrix scale-regions -S ../ChIP/H2Aub.bw  ../ChIP/H3K27ac.bw  ../ChIP/H3K27me3.bw  ../ChIP/H3K4me1.bw  ../ChIP/H3K4me3.bw  ../ChIP/H3K79me3.bw  ../ChIP/H3K9me3.bw \
-R MandB_ON.bed MandB_OFF.bed RD_MEF_IVF_NT_B.t.bed RU_MEF_IVF_NT_B.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName MB --sortRegions no

plotHeatmap -m MB \
-out MB_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' 

plotProfile -m MB \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' \
-out MB_profile.pdf

