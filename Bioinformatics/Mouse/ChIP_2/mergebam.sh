#cp SRR3208752_trimmed_unique.bam Input_2C.bam
samtools merge H3K4me3_2C.bam SRR3208752_trimmed_unique.bam SRR3208753_trimmed_unique.bam SRR3208754_trimmed_unique.bam -f
samtools merge H3K27me3_2C.bam SRR3208756_trimmed_unique.bam SRR3208757_trimmed_unique.bam SRR3208758_trimmed_unique.bam -f
#cp SRR3208759_trimmed_unique.bam Input_4C.bam
samtools merge H3K4me3_4C.bam SRR3208760_trimmed_unique.bam SRR3208761_trimmed_unique.bam SRR3208762_trimmed_unique.bam -f
samtools merge H3K27me3_4C.bam SRR3208763_trimmed_unique.bam SRR3208764_trimmed_unique.bam SRR3208765_trimmed_unique.bam -f
#cp SRR3208766_trimmed_unique.bam Input_8C.bam
samtools merge H3K4me3_8C.bam SRR3208767_trimmed_unique.bam SRR3208768_trimmed_unique.bam -f
samtools merge H3K27me3_8C.bam SRR3208769_trimmed_unique.bam SRR3208770_trimmed_unique.bam SRR3208771_trimmed_unique.bam -f
#cp SRR3208772_trimmed_unique.bam Input_M.bam
samtools merge H3K4me3_M.bam SRR3208773_trimmed_unique.bam SRR3208774_trimmed_unique.bam SRR3208775_trimmed_unique.bam -f
samtools merge H3K27me3_M.bam SRR3208776_trimmed_unique.bam SRR3208777_trimmed_unique.bam  -f
#cp SRR3208778_trimmed_unique.bam Input_B.bam
samtools merge H3K4me3_B.bam SRR3208779_trimmed_unique.bam SRR3208780_trimmed_unique.bam -f
samtools merge H3K27me3_B.bam SRR3208781_trimmed_unique.bam SRR3208782_trimmed_unique.bam -f

samtools index Input_2C.bam
samtools index Input_4C.bam
samtools index Input_8C.bam
samtools index Input_M.bam
samtools index Input_B.bam
samtools index H3K4me3_2C.bam
samtools index H3K4me3_4C.bam
samtools index H3K4me3_8C.bam
samtools index H3K4me3_M.bam
samtools index H3K4me3_B.bam
samtools index H3K27me3_2C.bam
samtools index H3K27me3_4C.bam
samtools index H3K27me3_8C.bam
samtools index H3K27me3_M.bam
samtools index H3K27me3_B.bam


bamCompare -b1 H3K4me3_2C.bam -b2 Input_2C.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K4me3_2C_vsInput.bw
bamCompare -b1 H3K4me3_4C.bam -b2 Input_4C.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K4me3_4C_vsInput.bw 
bamCompare -b1 H3K4me3_8C.bam -b2 Input_8C.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K4me3_8C_vsInput.bw 
bamCompare -b1 H3K4me3_M.bam -b2 Input_M.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K4me3_M_vsInput.bw 
bamCompare -b1 H3K4me3_B.bam -b2 Input_B.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K4me3_B_vsInput.bw 
bamCompare -b1 H3K27me3_2C.bam -b2 Input_2C.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K27me3_2C_vsInput.bw
bamCompare -b1 H3K27me3_4C.bam -b2 Input_4C.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K27me3_4C_vsInput.bw 
bamCompare -b1 H3K427me3_8C.bam -b2 Input_8C.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K27me3_8C_vsInput.bw 
bamCompare -b1 H3K27me3_M.bam -b2 Input_M.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K27me3_M_vsInput.bw 
bamCompare -b1 H3K27me3_B.bam -b2 Input_B.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K27me3_B_vsInput.bw 


#Peaklist : individual memory lists
#Conserved On and Off

computeMatrix scale-regions -S H3K4me3_2C_vsInput.bw H3K4me3_4C_vsInput.bw H3K4me3_8C_vsInput.bw H3K4me3_M_vsInput.bw H3K4me3_B_vsInput.bw \
-R PEAKSLIST
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName H3K4me3_plots --sortRegions no

computeMatrix scale-regions -S H3K27me3_2C_vsInput.bw H3K27me3_4C_vsInput.bw H3K27me3_8C_vsInput.bw H3K27me3_M_vsInput.bw H3K27me3_B_vsInput.bw \
-R PEAKSLIST
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 \
--binSize 10 --outFileName H3K27me3_plots --sortRegions no

plotHeatmap -m H3K4me3_plots  \
-out H3K4me3_plots.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Oranges Oranges Oranges \
--samplesLabel 2C 4C 8C Morula ICM \
--startLabel 'PS' \
--endLabel 'PE'


