#bamCompare -b1 H3K27me3_8C.bam -b2 Input_8C.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K27me3_8C_vsInput.bw 
#bamCompare -b1 H3K27me3_M.bam -b2 Input_M.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K27me3_M_vsInput.bw 
#bamCompare -b1 H3K27me3_B.bam -b2 Input_B.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K27me3_B_vsInput.bw 
#
#
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2/ SRR3208759.fastq
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR3208759_trimmed.fq -S SRR3208759_trimmed.sam
#samtools view -bS SRR3208759_trimmed.sam | samtools sort - SRR3208759_trimmed_sorted
#samtools index SRR3208759_trimmed_sorted.bam  SRR3208759_trimmed_sorted.bai
#samtools view -b -q 5 SRR3208759_trimmed_sorted.bam > SRR3208759_trimmed_unique.bam
#samtools index SRR3208759_trimmed_unique.bam SRR3208759_trimmed_unique.bai

#cp SRR3208759_trimmed_unique.bam Input_4C.bam

samtools index Input_4C.bam Input_4C.bai
bamCompare -b1 H3K4me3_4C.bam -b2 Input_4C.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K4me3_4C_vsInput.bw 
bamCompare -b1 H3K27me3_4C.bam -b2 Input_4C.bam --binSize 200 --extendReads 150 --normalizeUsing CPM -o H3K27me3_4C_vsInput.bw 

