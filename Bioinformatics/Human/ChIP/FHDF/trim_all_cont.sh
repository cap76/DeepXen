~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Human/ChIP/FHDF/ SRR227590.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Human/ChIP/FHDF/ SRR227591.fastq


bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/Bowtie2Index/genome -U SRR227591_trimmed.fq -S SRR227591_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/Bowtie2Index/genome -U SRR227590_trimmed.fq -S SRR227590_trimmed.sam

samtools view -bS	SRR227590_trimmed.sam	| samtools sort -	SRR227590_trimmed_sorted
samtools view -bS	SRR227591_trimmed.sam	| samtools sort -	SRR227591_trimmed_sorted
samtools index	SRR227590_trimmed_sorted.bam	SRR227590_trimmed_sorted.bai
samtools index	SRR227591_trimmed_sorted.bam	SRR227591_trimmed_sorted.bai
samtools view -b -q 5 	SRR227590_trimmed_sorted.bam >	SRR227590_trimmed_sorted_unique.bam
samtools view -b -q 5 	SRR227591_trimmed_sorted.bam >	SRR227591_trimmed_sorted_unique.bam
samtools index SRR227590_trimmed_sorted_unique.bam SRR227590_trimmed_sorted_unique.bai
samtools index SRR227591_trimmed_sorted_unique.bam SRR227591_trimmed_sorted_unique.bai
samtools merge input.bam SRR227590_trimmed_sorted.bam SRR227591_trimmed_sorted.bam
samtools index input.bam input.bai
samtools index H3K4me3.bam H3K4me3.bai
samtools index H3K9ac.bam H3K9ac.bai
samtools index H3K27ac.bam H3K27ac.bai
samtools index H3K27me3.bam H3K27me3.bai
samtools index H3K4me2.bam H3K4me2.bai
samtools index H4K20me1.bam H4K20me1.bai
samtools index H2A.Z.bam H2A.Z.bai
samtools index H3K4me1.bam H3K4me1.bai
samtools index H3K9me3.bam H3K9me3.bai

bamCoverage --bam H3K4me3.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K4me3.bw
bamCoverage --bam H3K9ac.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o  H3K9ac.bw
bamCoverage --bam H3K27ac.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K27ac.bw
bamCoverage --bam H3K36me3.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K36me3.bw
bamCoverage --bam H3K27me3.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K27me3.bw
bamCoverage --bam H3K4me2.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K4me2.bw
bamCoverage --bam H4K20me1.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o H4K20me1.bw
bamCoverage --bam H2A.Z.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o H2A.Z.bw
bamCoverage --bam H3K4me1.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K4me1.bw
bamCoverage --bam H3K9me3.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K9me3.bw


bamCompare --bam -b1 H3K4me3.bam -b2 input.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K4me3vInp.bw
bamCompare --bam -b1 H3K9ac.bam -b2 input.bam  --binSize 10  --extendReads 150  --normalizeUsing CPM -o  H3K9acvInp.bw
bamCompare --bam -b1 H3K27ac.bam -b2 input.bam  --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K27acvInp.bw
bamCompare --bam -b1 H3K36me3.bam -b2 input.bam  --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K36me3vInp.bw
bamCompare --bam -b1 H3K27me3.bam -b2 input.bam  --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K27me3vInp.bw
bamCompare --bam -b1 H3K4me2.bam -b2 input.bam  --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K4me2vInp.bw
bamCompare --bam -b1 H4K20me1.bam -b2 input.bam  --binSize 10  --extendReads 150  --normalizeUsing CPM -o H4K20me1vInp.bw
bamCompare --bam -b1 H2A.Z.bam -b2 input.bam  --binSize 10  --extendReads 150  --normalizeUsing CPM -o H2A.ZvInp.bw
bamCompare --bam -b1 H3K4me1.bam -b2 input.bam  --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K4me1vInp.bw
bamCompare --bam -b1 H3K9me3.bam -b2 input.bam  --binSize 10  --extendReads 150  --normalizeUsing CPM -o H3K9me3vInp.bw



#bamCoverage --bam SRR5029389_trimmed_sorted_unique.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o SRR5029389_trimmed_sorted_unique.bw 
#bamCoverage --bam SRR5029390_trimmed_sorted_unique.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o SRR5029390_trimmed_sorted_unique.bw 
#bamCoverage --bam SRR5029391_trimmed_sorted_unique.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o SRR5029391_trimmed_sorted_unique.bw 
#bamCoverage --bam SRR5029392_trimmed_sorted_unique.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o SRR5029392_trimmed_sorted_unique.bw 

#bamCoverage --bam 	SRR227372_trimmed_sorted_unique.bam	 --binSize 10  --extendReads 150  --normalizeUsing CPM -o 	SRR227372_trimmed_sorted_unique.bw
#bamCoverage --bam 	SRR227373_trimmed_sorted_unique.bam	 --binSize 10  --extendReads 150  --normalizeUsing CPM -o 	SRR227373_trimmed_sorted_unique.bw

