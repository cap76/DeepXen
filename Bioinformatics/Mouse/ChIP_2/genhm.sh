
computeMatrix scale-regions -S H3K27me3_2C_vsInput.bw H3K27me3_2C_vsInput.bw H3K27me3_8C_vsInput.bw H3K27me3_M_vsInput.bw H3K27me3_B_vsInput.bw H3K4me3_2C_vsInput.bw H3K4me3_4C_vsInput.bw H3K4me3_8C_vsInput.bw H3K4me3_M_vsInput.bw H3K4me3_B_vsInput.bw \
-R ../Mouse/RNAseq/ON_MEF_IVF_NT_1C.t.bed ../Mouse/RNAseq/OFF_MEF_IVF_NT_1C.t.bed ../Mouse/RNAseq/RD_MEF_IVF_NT_1C.t.bed ../Mouse/RNAseq/RU_MEF_IVF_NT_1C.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName 1C --sortRegions no

plotHeatmap -m 1C \
-out 1C_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--startLabel 'PS' \
--endLabel 'PE'

plotProfile -m 1C \
--startLabel 'PS' \
--endLabel 'PE' \
-out 1C_profile.pdf

computeMatrix scale-regions -S H3K27me3_2C_vsInput.bw H3K27me3_2C_vsInput.bw H3K27me3_8C_vsInput.bw H3K27me3_M_vsInput.bw H3K27me3_B_vsInput.bw H3K4me3_2C_vsInput.bw H3K4me3_4C_vsInput.bw H3K4me3_8C_vsInput.bw H3K4me3_M_vsInput.bw H3K4me3_B_vsInput.bw \
-R ../Mouse/RNAseq/ON_MEF_IVF_NT_2C.t.bed ../Mouse/RNAseq/OFF_MEF_IVF_NT_2C.t.bed ../Mouse/RNAseq/RD_MEF_IVF_NT_2C.t.bed ../Mouse/RNAseq/RU_MEF_IVF_NT_2C.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName 2C --sortRegions no

plotHeatmap -m 2C \
-out 2C_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--startLabel 'PS' \
--endLabel 'PE'

plotProfile -m 2C \
--startLabel 'PS' \
--endLabel 'PE' \
-out 2C_profile.pdf

computeMatrix scale-regions -S H3K27me3_2C_vsInput.bw H3K27me3_2C_vsInput.bw H3K27me3_8C_vsInput.bw H3K27me3_M_vsInput.bw H3K27me3_B_vsInput.bw H3K4me3_2C_vsInput.bw H3K4me3_4C_vsInput.bw H3K4me3_8C_vsInput.bw H3K4me3_M_vsInput.bw H3K4me3_B_vsInput.bw \
-R ../Mouse/RNAseq/ON_MEF_IVF_NT_4C.t.bed ../Mouse/RNAseq/OFF_MEF_IVF_NT_4C.t.bed ../Mouse/RNAseq/RD_MEF_IVF_NT_4C.t.bed ../Mouse/RNAseq/RU_MEF_IVF_NT_4C.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName 4C --sortRegions no

plotHeatmap -m 4C \
-out 4C_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--startLabel 'PS' \
--endLabel 'PE'

plotProfile -m 4C \
--startLabel 'PS' \
--endLabel 'PE' \
-out 4C_profile.pdf


computeMatrix scale-regions -S H3K27me3_2C_vsInput.bw H3K27me3_2C_vsInput.bw H3K27me3_8C_vsInput.bw H3K27me3_M_vsInput.bw H3K27me3_B_vsInput.bw H3K4me3_2C_vsInput.bw H3K4me3_4C_vsInput.bw H3K4me3_8C_vsInput.bw H3K4me3_M_vsInput.bw H3K4me3_B_vsInput.bw \
-R ../Mouse/RNAseq/ON_MEF_IVF_NT_8C.t.bed ../Mouse/RNAseq/OFF_MEF_IVF_NT_8C.t.bed ../Mouse/RNAseq/RD_MEF_IVF_NT_8C.t.bed ../Mouse/RNAseq/RU_MEF_IVF_NT_8C.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName 8C --sortRegions no

plotHeatmap -m 8C \
-out 8C_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--startLabel 'PS' \
--endLabel 'PE'

plotProfile -m 8C \
--startLabel 'PS' \
--endLabel 'PE' \
-out 8C_profile.pdf


computeMatrix scale-regions -S H3K27me3_2C_vsInput.bw H3K27me3_2C_vsInput.bw H3K27me3_8C_vsInput.bw H3K27me3_M_vsInput.bw H3K27me3_B_vsInput.bw H3K4me3_2C_vsInput.bw H3K4me3_4C_vsInput.bw H3K4me3_8C_vsInput.bw H3K4me3_M_vsInput.bw H3K4me3_B_vsInput.bw \
-R ../Mouse/RNAseq/ON_MEF_IVF_NT_B.t.bed ../Mouse/RNAseq/OFF_MEF_IVF_NT_B.t.bed ../Mouse/RNAseq/RD_MEF_IVF_NT_B.t.bed ../Mouse/RNAseq/RU_MEF_IVF_NT_B.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName BC --sortRegions no

plotHeatmap -m BC \
-out BC_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--startLabel 'PS' \
--endLabel 'PE'

plotProfile -m BC \
--startLabel 'PS' \
--endLabel 'PE' \
-out BC_profile.pdf


computeMatrix scale-regions -S H3K27me3_2C_vsInput.bw H3K27me3_2C_vsInput.bw H3K27me3_8C_vsInput.bw H3K27me3_M_vsInput.bw H3K27me3_B_vsInput.bw H3K4me3_2C_vsInput.bw H3K4me3_4C_vsInput.bw H3K4me3_8C_vsInput.bw H3K4me3_M_vsInput.bw H3K4me3_B_vsInput.bw \
-R ../Mouse/RNAseq/ON_MEF_IVF_NT_M.t.bed ../Mouse/RNAseq/OFF_MEF_IVF_NT_M.t.bed ../Mouse/RNAseq/RD_MEF_IVF_NT_M.t.bed ../Mouse/RNAseq/RU_MEF_IVF_NT_M.t.bed \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName MC --sortRegions no

plotHeatmap -m MC \
-out MC_heatmap.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Reds Greens Blues Oranges Reds Greens Blues \
--startLabel 'PS' \
--endLabel 'PE'

plotProfile -m MC \
--startLabel 'PS' \
--endLabel 'PE' \
-out MC_profile.pdf

