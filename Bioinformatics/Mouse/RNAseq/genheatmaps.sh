computeMatrix scale-regions -S ../ChIP/H2Aub.bw  ../ChIP/H3K27ac.bw  ../ChIP/H3K27me3.bw  ../ChIP/H3K4me1.bw  ../ChIP/H3K4me3.bw  ../ChIP/H3K79me3.bw  ../ChIP/H3K9me3.bw \
-R ON_MEF_IVF_NT_1C.t.bed OFF_MEF_IVF_NT_1C.t.bed RD_MEF_IVF_NT_1C.t.bed RU_MEF_IVF_NT_1C.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName 1C --sortRegions no

plotHeatmap -m 1C \
-out 1C_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE'

plotProfile -m 1C \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' \
-out 1C_profile.pdf

computeMatrix scale-regions -S ../ChIP/H2Aub.bw  ../ChIP/H3K27ac.bw  ../ChIP/H3K27me3.bw  ../ChIP/H3K4me1.bw  ../ChIP/H3K4me3.bw  ../ChIP/H3K79me3.bw  ../ChIP/H3K9me3.bw \
-R ON_MEF_IVF_NT_2C.t.bed OFF_MEF_IVF_NT_2C.t.bed RD_MEF_IVF_NT_2C.t.bed RU_MEF_IVF_NT_2C.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName 2C --sortRegions no

plotHeatmap -m 2C \
-out 2C_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' 

plotProfile -m 2C \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' \
-out 2C_profile.pdf

computeMatrix scale-regions -S ../ChIP/H2Aub.bw  ../ChIP/H3K27ac.bw  ../ChIP/H3K27me3.bw  ../ChIP/H3K4me1.bw  ../ChIP/H3K4me3.bw  ../ChIP/H3K79me3.bw  ../ChIP/H3K9me3.bw \
-R ON_MEF_IVF_NT_4C.t.bed OFF_MEF_IVF_NT_4C.t.bed RD_MEF_IVF_NT_4C.t.bed RU_MEF_IVF_NT_4C.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName 4C --sortRegions no

plotHeatmap -m 4C \
-out 4C_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' 

plotProfile -m 4C \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' \
-out 4C_profile.pdf


computeMatrix scale-regions -S ../ChIP/H2Aub.bw  ../ChIP/H3K27ac.bw  ../ChIP/H3K27me3.bw  ../ChIP/H3K4me1.bw  ../ChIP/H3K4me3.bw  ../ChIP/H3K79me3.bw  ../ChIP/H3K9me3.bw \
-R ON_MEF_IVF_NT_8C.t.bed OFF_MEF_IVF_NT_8C.t.bed RD_MEF_IVF_NT_8C.t.bed RU_MEF_IVF_NT_8C.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName 8C --sortRegions no

plotHeatmap -m 8C \
-out 8C_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' 

plotProfile -m 8C \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' \
-out 8C_profile.pdf


computeMatrix scale-regions -S ../ChIP/H2Aub.bw  ../ChIP/H3K27ac.bw  ../ChIP/H3K27me3.bw  ../ChIP/H3K4me1.bw  ../ChIP/H3K4me3.bw  ../ChIP/H3K79me3.bw  ../ChIP/H3K9me3.bw \
-R ON_MEF_IVF_NT_M.t.bed OFF_MEF_IVF_NT_M.t.bed RD_MEF_IVF_NT_M.t.bed RU_MEF_IVF_NT_M.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName MC --sortRegions no

plotHeatmap -m MC \
-out MC_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' 

plotProfile -m MC \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' \
-out MC_profile.pdf

computeMatrix scale-regions -S ../ChIP/H2Aub.bw  ../ChIP/H3K27ac.bw  ../ChIP/H3K27me3.bw  ../ChIP/H3K4me1.bw  ../ChIP/H3K4me3.bw  ../ChIP/H3K79me3.bw  ../ChIP/H3K9me3.bw \
-R ON_MEF_IVF_NT_B.t.bed OFF_MEF_IVF_NT_B.t.bed RD_MEF_IVF_NT_B.t.bed RU_MEF_IVF_NT_B.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName BC --sortRegions no

plotHeatmap -m 1C \
-out BC_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' 

plotProfile -m BC \
--samplesLabel H2Aub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K79me3.bw H3K9me3.bw \
--startLabel 'PS' \
--endLabel 'PE' \
-out BC_profile.pdf




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
