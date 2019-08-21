#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Human/ChIP SRR5029389.fastq  
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Human/ChIP SRR5029390.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Human/ChIP SRR5029391.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Human/ChIP SRR5029392.fastq

#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/Bowtie2Index/genome -U SRR5029389_trimmed.fq -S SRR5029389_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/Bowtie2Index/genome -U SRR5029390_trimmed.fq -S SRR5029390_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/Bowtie2Index/genome -U SRR5029391_trimmed.fq -S SRR5029391_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Sequence/Bowtie2Index/genome -U SRR5029392_trimmed.fq -S SRR5029392_trimmed.sam


samtools view -bS SRR5029389_trimmed.sam | samtools sort - SRR5029389_trimmed_sorted
samtools view -bS SRR5029390_trimmed.sam | samtools sort - SRR5029390_trimmed_sorted
samtools view -bS SRR5029391_trimmed.sam | samtools sort - SRR5029391_trimmed_sorted
samtools view -bS SRR5029392_trimmed.sam | samtools sort - SRR5029392_trimmed_sorted

samtools index SRR5029389_trimmed_sorted.bam SRR5029389_trimmed_sorted.bai
samtools index SRR5029390_trimmed_sorted.bam SRR5029390_trimmed_sorted.bai
samtools index SRR5029391_trimmed_sorted.bam SRR5029391_trimmed_sorted.bai
samtools index SRR5029392_trimmed_sorted.bam SRR5029392_trimmed_sorted.bai


samtools view -b -q 5 SRR5029389_trimmed_sorted.bam > SRR5029389_trimmed_sorted_unique.bam
samtools view -b -q 5 SRR5029390_trimmed_sorted.bam > SRR5029390_trimmed_sorted_unique.bam
samtools view -b -q 5 SRR5029391_trimmed_sorted.bam > SRR5029391_trimmed_sorted_unique.bam
samtools view -b -q 5 SRR5029392_trimmed_sorted.bam > SRR5029392_trimmed_sorted_unique.bam

samtools index SRR5029389_trimmed_sorted_unique.bam SRR5029389_trimmed_sorted_unique.bai
samtools index SRR5029390_trimmed_sorted_unique.bam SRR5029390_trimmed_sorted_unique.bai
samtools index SRR5029391_trimmed_sorted_unique.bam SRR5029391_trimmed_sorted_unique.bai
samtools index SRR5029392_trimmed_sorted_unique.bam SRR5029392_trimmed_sorted_unique.bai

bamCoverage --bam SRR5029389_trimmed_sorted_unique.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o SRR5029389_trimmed_sorted_unique.bw 
bamCoverage --bam SRR5029390_trimmed_sorted_unique.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o SRR5029390_trimmed_sorted_unique.bw 
bamCoverage --bam SRR5029391_trimmed_sorted_unique.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o SRR5029391_trimmed_sorted_unique.bw 
bamCoverage --bam SRR5029392_trimmed_sorted_unique.bam --binSize 10  --extendReads 150  --normalizeUsing CPM -o SRR5029392_trimmed_sorted_unique.bw 

