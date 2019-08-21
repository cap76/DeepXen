
#2 cell input
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/002/SRR3208752/SRR3208752_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/002/SRR3208752/SRR3208752_2.fastq.gz
#H3K4
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/003/SRR3208753/SRR3208753_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/003/SRR3208753/SRR3208753_2.fastq.gz

#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/004/SRR3208754/SRR3208754_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/004/SRR3208754/SRR3208754_2.fastq.gz

#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/005/SRR3208755/SRR3208755_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/005/SRR3208755/SRR3208755_2.fastq.gz

#H3K27me3
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/006/SRR3208756/SRR3208756_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/006/SRR3208756/SRR3208756_2.fastq.gz

#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/007/SRR3208757/SRR3208757_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/007/SRR3208757/SRR3208757_2.fastq.gz

#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/008/SRR3208758/SRR3208758_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/008/SRR3208758/SRR3208758_2.fastq.gz


#cp SRR3208752_trimmed_unique.bam Input_2C.bam
#samtools merge H3K4me3_2C.bam SRR3208752_trimmed_unique.bam SRR3208753_trimmed_unique.bam SRR3208754_trimmed_unique.bam
#samtools merge H3K27me3_2C.bam SRR3208756_trimmed_unique.bam SRR3208757_trimmed_unique.bam SRR3208758_trimmed_unique.bam

~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208752_1.fastq SRR3208752_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208753_1.fastq SRR3208753_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208754_1.fastq SRR3208754_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208755_1.fastq SRR3208755_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208756_1.fastq SRR3208756_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208757_1.fastq SRR3208757_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208758_1.fastq SRR3208758_2.fastq

bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208752_1_val_1.fq -2 SRR3208752_2_val_2.fq -S SRR3208752_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208753_1_val_1.fq -2 SRR3208753_2_val_2.fq -S SRR3208753_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208754_1_val_1.fq -2 SRR3208754_2_val_2.fq -S SRR3208754_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208755_1_val_1.fq -2 SRR3208755_2_val_2.fq -S SRR3208755_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208756_1_val_1.fq -2 SRR3208756_2_val_2.fq -S SRR3208756_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208757_1_val_1.fq -2 SRR3208757_2_val_2.fq -S SRR3208757_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208758_1_val_1.fq -2 SRR3208758_2_val_2.fq -S SRR3208758_trimmed.sam

samtools view -bS SRR3208752_trimmed.sam | samtools sort - SRR3208752_trimmed_sorted
samtools index SRR3208752_trimmed_sorted.bam  SRR3208752_trimmed_sorted.bai
samtools view -b -q 5 SRR3208752_trimmed_sorted.bam > SRR3208752_trimmed_unique.bam
samtools index SRR3208752_trimmed_unique.bam SRR3208752_trimmed_unique.bai

samtools view -bS SRR3208753_trimmed.sam | samtools sort - SRR3208753_trimmed_sorted
samtools index SRR3208753_trimmed_sorted.bam  SRR3208753_trimmed_sorted.bai
samtools view -b -q 5 SRR3208753_trimmed_sorted.bam > SRR3208753_trimmed_unique.bam
samtools index SRR3208753_trimmed_unique.bam SRR3208753_trimmed_unique.bai

samtools view -bS SRR3208754_trimmed.sam | samtools sort - SRR3208754_trimmed_sorted
samtools index SRR3208754_trimmed_sorted.bam  SRR3208754_trimmed_sorted.bai
samtools view -b -q 5 SRR3208754_trimmed_sorted.bam > SRR3208754_trimmed_unique.bam
samtools index SRR3208754_trimmed_unique.bam SRR3208754_trimmed_unique.bai

samtools view -bS SRR3208755_trimmed.sam | samtools sort - SRR3208755_trimmed_sorted
samtools index SRR3208755_trimmed_sorted.bam  SRR3208755_trimmed_sorted.bai
samtools view -b -q 5 SRR32087525_trimmed_sorted.bam > SRR3208755_trimmed_unique.bam
samtools index SRR3208755_trimmed_unique.bam SRR3208755_trimmed_unique.bai

samtools view -bS SRR3208756_trimmed.sam | samtools sort - SRR3208756_trimmed_sorted
samtools index SRR3208756_trimmed_sorted.bam  SRR3208756_trimmed_sorted.bai
samtools view -b -q 5 SRR3208756_trimmed_sorted.bam > SRR3208756_trimmed_unique.bam
samtools index SRR3208756_trimmed_unique.bam SRR3208756_trimmed_unique.bai

samtools view -bS SRR3208757_trimmed.sam | samtools sort - SRR3208757_trimmed_sorted
samtools index SRR3208757_trimmed_sorted.bam  SRR3208757_trimmed_sorted.bai
samtools view -b -q 5 SRR3208757_trimmed_sorted.bam > SRR3208757_trimmed_unique.bam
samtools index SRR3208757_trimmed_unique.bam SRR3208757_trimmed_unique.bai

samtools view -bS SRR3208758_trimmed.sam | samtools sort - SRR3208758_trimmed_sorted
samtools index SRR3208758_trimmed_sorted.bam  SRR3208758_trimmed_sorted.bai
samtools view -b -q 5 SRR3208758_trimmed_sorted.bam > SRR3208758_trimmed_unique.bam
samtools index SRR3208758_trimmed_unique.bam SRR3208758_trimmed_unique.bai


#4cell
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/009/SRR3208759/SRR3208759.fastq.gz
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2/ SRR3208759..fastq
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR3208759_trimmed.fq -S SRR3208759_trimmed.sam
samtools view -bS SRR3208759_trimmed.sam | samtools sort - SRR3208759_trimmed_sorted
samtools index SRR3208759_trimmed_sorted.bam  SRR3208759_trimmed_sorted.bai
samtools view -b -q 5 SRR3208759_trimmed_sorted.bam > SRR3208759_trimmed_unique.bam
samtools index SRR3208759_trimmed_unique.bam SRR3208759_trimmed_unique.bai



#H3K4me3
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/000/SRR3208760/SRR3208760_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/001/SRR3208761/SRR3208761_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/002/SRR3208762/SRR3208762_1.fastq.gz
#H3K27me3
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/003/SRR3208763/SRR3208763_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/004/SRR3208764/SRR3208764_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/005/SRR3208765/SRR3208765_1.fastq.gz

#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/000/SRR3208760/SRR3208760_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/001/SRR3208761/SRR3208761_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/002/SRR3208762/SRR3208762_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/003/SRR3208763/SRR3208763_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/004/SRR3208764/SRR3208764_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/005/SRR3208765/SRR3208765_2.fastq.gz


#cp SRR3208752_trimmed_unique.bam Input_2C.bam
#samtools merge H3K4me3_2C.bam SRR3208752_trimmed_unique.bam SRR3208753_trimmed_unique.bam SRR3208754_trimmed_unique.bam
#samtools merge H3K27me3_2C.bam SRR3208756_trimmed_unique.bam SRR3208757_trimmed_unique.bam SRR3208758_trimmed_unique.bam
#cp SRR3208759_trimmed_unique.bam Input_4C.bam
#samtools merge H3K4me3_4C.bam SRR3208760_trimmed_unique.bam SRR3208761_trimmed_unique.bam SRR3208762_trimmed_unique.bam
#samtools merge H3K27me3_4C.bam SRR3208763_trimmed_unique.bam SRR3208764_trimmed_unique.bam SRR3208765_trimmed_unique.bam



~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208760_1.fastq SRR3208760_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208761_1.fastq SRR3208761_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208762_1.fastq SRR3208762_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208763_1.fastq SRR3208763_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208764_1.fastq SRR3208764_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208765_1.fastq SRR3208765_2.fastq


bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208760_1_val_1.fq -2 SRR3208760_2_val_2.fq -S SRR3208760_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208761_1_val_1.fq -2 SRR3208761_2_val_2.fq -S SRR3208761_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208762_1_val_1.fq -2 SRR3208762_2_val_2.fq -S SRR3208762_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208763_1_val_1.fq -2 SRR3208763_2_val_2.fq -S SRR3208763_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208764_1_val_1.fq -2 SRR3208764_2_val_2.fq -S SRR3208764_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208765_1_val_1.fq -2 SRR3208765_2_val_2.fq -S SRR3208765_trimmed.sam



samtools view -bS SRR3208760_trimmed.sam | samtools sort - SRR3208760_trimmed_sorted
samtools index SRR3208760_trimmed_sorted.bam  SRR3208760_trimmed_sorted.bai
samtools view -b -q 5 SRR3208760_trimmed_sorted.bam > SRR3208760_trimmed_unique.bam
samtools index SRR3208760_trimmed_unique.bam SRR3208760_trimmed_unique.bai

samtools view -bS SRR3208761_trimmed.sam | samtools sort - SRR3208761_trimmed_sorted
samtools index SRR3208761_trimmed_sorted.bam  SRR3208761_trimmed_sorted.bai
samtools view -b -q 5 SRR3208761_trimmed_sorted.bam > SRR3208761_trimmed_unique.bam
samtools index SRR3208761_trimmed_unique.bam SRR3208761_trimmed_unique.bai

samtools view -bS SRR3208762_trimmed.sam | samtools sort - SRR3208762_trimmed_sorted
samtools index SRR3208762_trimmed_sorted.bam  SRR3208762_trimmed_sorted.bai
samtools view -b -q 5 SRR3208762_trimmed_sorted.bam > SRR3208762_trimmed_unique.bam
samtools index SRR3208762_trimmed_unique.bam SRR3208762_trimmed_unique.bai

samtools view -bS SRR3208763_trimmed.sam | samtools sort - SRR3208763_trimmed_sorted
samtools index SRR3208763_trimmed_sorted.bam  SRR3208763_trimmed_sorted.bai
samtools view -b -q 5 SRR3208763_trimmed_sorted.bam > SRR3208763_trimmed_unique.bam
samtools index SRR3208763_trimmed_unique.bam SRR3208763_trimmed_unique.bai

samtools view -bS SRR3208764_trimmed.sam | samtools sort - SRR3208764_trimmed_sorted
samtools index SRR3208764_trimmed_sorted.bam  SRR3208764_trimmed_sorted.bai
samtools view -b -q 5 SRR3208764_trimmed_sorted.bam > SRR3208764_trimmed_unique.bam
samtools index SRR3208764_trimmed_unique.bam SRR3208764_trimmed_unique.bai


samtools view -bS SRR3208765_trimmed.sam | samtools sort - SRR3208765_trimmed_sorted
samtools index SRR3208765_trimmed_sorted.bam  SRR3208765_trimmed_sorted.bai
samtools view -b -q 5 SRR3208765_trimmed_sorted.bam > SRR3208765_trimmed_unique.bam
samtools index SRR3208765_trimmed_unique.bam SRR3208765_trimmed_unique.bai




#8 cell
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/006/SRR3208766/SRR3208766_1.fastq.gz

#H3K4
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/007/SRR3208767/SRR3208767_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/008/SRR3208768/SRR3208768_1.fastq.gz
#K3K27me
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/009/SRR3208769/SRR3208769_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/000/SRR3208770/SRR3208770_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/001/SRR3208771/SRR3208771_1.fastq.gz

#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/006/SRR3208766/SRR3208766_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/007/SRR3208767/SRR3208767_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/008/SRR3208768/SRR3208768_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/009/SRR3208769/SRR3208769_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/000/SRR3208770/SRR3208770_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/001/SRR3208771/SRR3208771_2.fastq.gz


#cp SRR3208752_trimmed_unique.bam Input_2C.bam
#samtools merge H3K4me3_2C.bam SRR3208752_trimmed_unique.bam SRR3208753_trimmed_unique.bam SRR3208754_trimmed_unique.bam
#samtools merge H3K27me3_2C.bam SRR3208756_trimmed_unique.bam SRR3208757_trimmed_unique.bam SRR3208758_trimmed_unique.bam
#cp SRR3208759_trimmed_unique.bam Input_4C.bam
#samtools merge H3K4me3_4C.bam SRR3208760_trimmed_unique.bam SRR3208761_trimmed_unique.bam SRR3208762_trimmed_unique.bam
#samtools merge H3K27me3_4C.bam SRR3208763_trimmed_unique.bam SRR3208764_trimmed_unique.bam SRR3208765_trimmed_unique.bam
#cp SRR3208766_trimmed_unique.bam Input_8C.bam
#samtools merge H3K4me3_8C.bam SRR3208767_trimmed_unique.bam SRR3208768_trimmed_unique.bam 
#samtools merge H3K27me3_8C.bam SRR3208769_trimmed_unique.bam SRR3208770_trimmed_unique.bam SRR3208771_trimmed_unique.bam



~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208766_1.fastq SRR3208766_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208767_1.fastq SRR3208767_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208768_1.fastq SRR3208768_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208769_1.fastq SRR3208769_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208770_1.fastq SRR3208770_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208771_1.fastq SRR3208771_2.fastq


bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208766_1_val_1.fq -2 SRR3208766_2_val_2.fq -S SRR3208766_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208767_1_val_1.fq -2 SRR3208767_2_val_2.fq -S SRR3208767_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208768_1_val_1.fq -2 SRR3208768_2_val_2.fq -S SRR3208768_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208769_1_val_1.fq -2 SRR3208769_2_val_2.fq -S SRR3208769_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208770_1_val_1.fq -2 SRR3208770_2_val_2.fq -S SRR3208770_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208771_1_val_1.fq -2 SRR3208771_2_val_2.fq -S SRR3208771_trimmed.sam





samtools view -bS SRR3208766_trimmed.sam | samtools sort - SRR3208766_trimmed_sorted
samtools index SRR3208766_trimmed_sorted.bam  SRR3208766_trimmed_sorted.bai
samtools view -b -q 5 SRR3208766_trimmed_sorted.bam > SRR3208766_trimmed_unique.bam
samtools index SRR3208766_trimmed_unique.bam SRR3208766_trimmed_unique.bai

samtools view -bS SRR3208767_trimmed.sam | samtools sort - SRR3208767_trimmed_sorted
samtools index SRR3208767_trimmed_sorted.bam  SRR3208767_trimmed_sorted.bai
samtools view -b -q 5 SRR3208767_trimmed_sorted.bam > SRR3208767_trimmed_unique.bam
samtools index SRR3208767_trimmed_unique.bam SRR3208767_trimmed_unique.bai

samtools view -bS SRR3208768_trimmed.sam | samtools sort - SRR3208768_trimmed_sorted
samtools index SRR3208768_trimmed_sorted.bam  SRR3208768_trimmed_sorted.bai
samtools view -b -q 5 SRR3208768_trimmed_sorted.bam > SRR3208768_trimmed_unique.bam
samtools index SRR3208768_trimmed_unique.bam SRR3208768_trimmed_unique.bai

samtools view -bS SRR3208769_trimmed.sam | samtools sort - SRR3208769_trimmed_sorted
samtools index SRR3208769_trimmed_sorted.bam  SRR3208769_trimmed_sorted.bai
samtools view -b -q 5 SRR3208769_trimmed_sorted.bam > SRR3208769_trimmed_unique.bam
samtools index SRR3208769_trimmed_unique.bam SRR3208769_trimmed_unique.bai

samtools view -bS SRR3208770_trimmed.sam | samtools sort - SRR3208770_trimmed_sorted
samtools index SRR3208770_trimmed_sorted.bam  SRR3208770_trimmed_sorted.bai
samtools view -b -q 5 SRR3208770_trimmed_sorted.bam > SRR3208770_trimmed_unique.bam
samtools index SRR3208770_trimmed_unique.bam SRR3208770_trimmed_unique.bai

samtools view -bS SRR3208771_trimmed.sam | samtools sort - SRR3208771_trimmed_sorted
samtools index SRR3208771_trimmed_sorted.bam  SRR3208771_trimmed_sorted.bai
samtools view -b -q 5 SRR3208771_trimmed_sorted.bam > SRR3208771_trimmed_unique.bam
samtools index SRR3208771_trimmed_unique.bam SRR3208771_trimmed_unique.bai



#Morula
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/002/SRR3208772/SRR3208772_1.fastq.gz
#H3K4
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/003/SRR3208773/SRR3208773_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/004/SRR3208774/SRR3208774_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/005/SRR3208775/SRR3208775_1.fastq.gz
#H3K27
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/006/SRR3208776/SRR3208776_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/007/SRR3208777/SRR3208777_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/008/SRR3208778/SRR3208778_1.fastq.gz

#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/002/SRR3208772/SRR3208772_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/003/SRR3208773/SRR3208773_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/004/SRR3208774/SRR3208774_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/005/SRR3208775/SRR3208775_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/006/SRR3208776/SRR3208776_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/007/SRR3208777/SRR3208777_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/008/SRR3208778/SRR3208778_2.fastq.gz

#cp SRR3208752_trimmed_unique.bam Input_2C.bam
#samtools merge H3K4me3_2C.bam SRR3208752_trimmed_unique.bam SRR3208753_trimmed_unique.bam SRR3208754_trimmed_unique.bam
#samtools merge H3K27me3_2C.bam SRR3208756_trimmed_unique.bam SRR3208757_trimmed_unique.bam SRR3208758_trimmed_unique.bam
#cp SRR3208759_trimmed_unique.bam Input_4C.bam
#samtools merge H3K4me3_4C.bam SRR3208760_trimmed_unique.bam SRR3208761_trimmed_unique.bam SRR3208762_trimmed_unique.bam
#samtools merge H3K27me3_4C.bam SRR3208763_trimmed_unique.bam SRR3208764_trimmed_unique.bam SRR3208765_trimmed_unique.bam
#cp SRR3208766_trimmed_unique.bam Input_8C.bam
#samtools merge H3K4me3_8C.bam SRR3208767_trimmed_unique.bam SRR3208768_trimmed_unique.bam 
#samtools merge H3K27me3_8C.bam SRR3208769_trimmed_unique.bam SRR3208770_trimmed_unique.bam SRR3208771_trimmed_unique.bam
#cp SRR3208772_trimmed_unique.bam Input_M.bam
#samtools merge H3K4me3_M.bam SRR3208773_trimmed_unique.bam SRR3208774_trimmed_unique.bam SRR3208775_trimmed_unique.bam
#samtools merge H3K27me3_M.bam SRR3208776_trimmed_unique.bam SRR3208777_trimmed_unique.bam SRR3208778_trimmed_unique.bam



~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208772_1.fastq SRR3208772_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208773_1.fastq SRR3208773_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208774_1.fastq SRR3208774_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208775_1.fastq SRR3208775_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208776_1.fastq SRR3208776_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208777_1.fastq SRR3208777_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208778_1.fastq SRR3208778_2.fastq


bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208772_1_val_1.fq -2 SRR3208772_2_val_2.fq -S SRR3208772_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208773_1_val_1.fq -2 SRR3208773_2_val_2.fq -S SRR3208773_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208774_1_val_1.fq -2 SRR3208774_2_val_2.fq -S SRR3208774_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208775_1_val_1.fq -2 SRR3208775_2_val_2.fq -S SRR3208775_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208776_1_val_1.fq -2 SRR3208776_2_val_2.fq -S SRR3208776_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208777_1_val_1.fq -2 SRR3208777_2_val_2.fq -S SRR3208777_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208778_1_val_1.fq -2 SRR3208778_2_val_2.fq -S SRR3208778_trimmed.sam

samtools view -bS SRR3208772_trimmed.sam | samtools sort - SRR3208772_trimmed_sorted
samtools index SRR3208772_trimmed_sorted.bam  SRR3208772_trimmed_sorted.bai
samtools view -b -q 5 SRR3208772_trimmed_sorted.bam > SRR3208772_trimmed_unique.bam
samtools index SRR3208772_trimmed_unique.bam SRR3208772_trimmed_unique.bai

samtools view -bS SRR3208773_trimmed.sam | samtools sort - SRR3208773_trimmed_sorted
samtools index SRR3208773_trimmed_sorted.bam  SRR3208773_trimmed_sorted.bai
samtools view -b -q 5 SRR3208773_trimmed_sorted.bam > SRR3208773_trimmed_unique.bam
samtools index SRR3208773_trimmed_unique.bam SRR3208773_trimmed_unique.bai

samtools view -bS SRR3208774_trimmed.sam | samtools sort - SRR3208774_trimmed_sorted
samtools index SRR3208774_trimmed_sorted.bam  SRR3208774_trimmed_sorted.bai
samtools view -b -q 5 SRR3208774_trimmed_sorted.bam > SRR3208774_trimmed_unique.bam
samtools index SRR3208774_trimmed_unique.bam SRR3208774_trimmed_unique.bai

samtools view -bS SRR3208775_trimmed.sam | samtools sort - SRR3208775_trimmed_sorted
samtools index SRR3208775_trimmed_sorted.bam  SRR3208775_trimmed_sorted.bai
samtools view -b -q 5 SRR3208775_trimmed_sorted.bam > SRR3208775_trimmed_unique.bam
samtools index SRR3208775_trimmed_unique.bam SRR3208775_trimmed_unique.bai

samtools view -bS SRR3208776_trimmed.sam | samtools sort - SRR3208776_trimmed_sorted
samtools index SRR3208776_trimmed_sorted.bam  SRR3208776_trimmed_sorted.bai
samtools view -b -q 5 SRR3208776_trimmed_sorted.bam > SRR3208776_trimmed_unique.bam
samtools index SRR3208776_trimmed_unique.bam SRR3208776_trimmed_unique.bai

samtools view -bS SRR3208777_trimmed.sam | samtools sort - SRR3208777_trimmed_sorted
samtools index SRR3208777_trimmed_sorted.bam  SRR3208777_trimmed_sorted.bai
samtools view -b -q 5 SRR3208777_trimmed_sorted.bam > SRR3208777_trimmed_unique.bam
samtools index SRR3208777_trimmed_unique.bam SRR3208777_trimmed_unique.bai

samtools view -bS SRR3208778_trimmed.sam | samtools sort - SRR3208778_trimmed_sorted
samtools index SRR3208778_trimmed_sorted.bam  SRR3208778_trimmed_sorted.bai
samtools view -b -q 5 SRR3208778_trimmed_sorted.bam > SRR3208778_trimmed_unique.bam
samtools index SRR3208778_trimmed_unique.bam SRR3208778_trimmed_unique.bai

#ICM
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/008/SRR3208778/SRR3208778_1.fastq.gz
#H3K4me3
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/009/SRR3208779/SRR3208779_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/000/SRR3208780/SRR3208780_1.fastq.gz
#H3K27me3
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/001/SRR3208781/SRR3208781_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/002/SRR3208782/SRR3208782_1.fastq.gz


#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/008/SRR3208778/SRR3208778_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/009/SRR3208779/SRR3208779_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/000/SRR3208780/SRR3208780_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/001/SRR3208781/SRR3208781_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/002/SRR3208782/SRR3208782_2.fastq.gz


#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208778_1.fastq SRR3208778_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208779_1.fastq SRR3208779_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208780_1.fastq SRR3208780_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208781_1.fastq SRR3208781_2.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse_2 --paired SRR3208782_1.fastq SRR3208782_2.fastq


#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208778_1_val_1.fq -2 SRR3208778_2_val_2.fq -S SRR3208778_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208779_1_val_1.fq -2 SRR3208779_2_val_2.fq -S SRR3208779_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208780_1_val_1.fq -2 SRR3208780_2_val_2.fq -S SRR3208780_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208781_1_val_1.fq -2 SRR3208781_2_val_2.fq -S SRR3208781_trimmed.sam
bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -1 SRR3208782_1_val_1.fq -2 SRR3208782_2_val_2.fq -S SRR3208782_trimmed.sam


samtools view -bS SRR3208779_trimmed.sam | samtools sort - SRR3208779_trimmed_sorted
samtools index SRR3208779_trimmed_sorted.bam  SRR3208779_trimmed_sorted.bai
samtools view -b -q 5 SRR3208779_trimmed_sorted.bam > SRR3208779_trimmed_unique.bam
samtools index SRR3208779_trimmed_unique.bam SRR3208779_trimmed_unique.bai

samtools view -bS SRR3208780_trimmed.sam | samtools sort - SRR3208780_trimmed_sorted
samtools index SRR3208780_trimmed_sorted.bam  SRR3208780_trimmed_sorted.bai
samtools view -b -q 5 SRR3208780_trimmed_sorted.bam > SRR3208780_trimmed_unique.bam
samtools index SRR3208780_trimmed_unique.bam SRR3208780_trimmed_unique.bai

samtools view -bS SRR3208781_trimmed.sam | samtools sort - SRR3208781_trimmed_sorted
samtools index SRR3208781_trimmed_sorted.bam  SRR3208781_trimmed_sorted.bai
samtools view -b -q 5 SRR3208781_trimmed_sorted.bam > SRR3208781_trimmed_unique.bam
samtools index SRR3208781_trimmed_unique.bam SRR3208781_trimmed_unique.bai

samtools view -bS SRR3208782_trimmed.sam | samtools sort - SRR3208782_trimmed_sorted
samtools index SRR3208782_trimmed_sorted.bam  SRR3208782_trimmed_sorted.bai
samtools view -b -q 5 SRR3208782_trimmed_sorted.bam > SRR3208782_trimmed_unique.bam
samtools index SRR3208782_trimmed_unique.bam SRR3208782_trimmed_unique.bai


cp SRR3208752_trimmed_unique.bam Input_2C.bam
samtools merge H3K4me3_2C.bam SRR3208752_trimmed_unique.bam SRR3208753_trimmed_unique.bam SRR3208754_trimmed_unique.bam
samtools merge H3K27me3_2C.bam SRR3208756_trimmed_unique.bam SRR3208757_trimmed_unique.bam SRR3208758_trimmed_unique.bam
cp SRR3208759_trimmed_unique.bam Input_4C.bam
samtools merge H3K4me3_4C.bam SRR3208760_trimmed_unique.bam SRR3208761_trimmed_unique.bam SRR3208762_trimmed_unique.bam
samtools merge H3K27me3_4C.bam SRR3208763_trimmed_unique.bam SRR3208764_trimmed_unique.bam SRR3208765_trimmed_unique.bam
cp SRR3208766_trimmed_unique.bam Input_8C.bam
samtools merge H3K4me3_8C.bam SRR3208767_trimmed_unique.bam SRR3208768_trimmed_unique.bam 
samtools merge H3K27me3_8C.bam SRR3208769_trimmed_unique.bam SRR3208770_trimmed_unique.bam SRR3208771_trimmed_unique.bam
cp SRR3208772_trimmed_unique.bam Input_M.bam
samtools merge H3K4me3_M.bam SRR3208773_trimmed_unique.bam SRR3208774_trimmed_unique.bam SRR3208775_trimmed_unique.bam
samtools merge H3K27me3_M.bam SRR3208776_trimmed_unique.bam SRR3208777_trimmed_unique.bam 
cp SRR3208778_trimmed_unique.bam Input_B.bam
samtools merge H3K4me3_B.bam SRR3208779_trimmed_unique.bam SRR3208780_trimmed_unique.bam 
samtools merge H3K27me3_B.bam SRR3208781_trimmed_unique.bam SRR3208782_trimmed_unique.bam

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

#--samplesLabel Brd4_ES Brd4_EpiLC Epi_1 Epi_2 Epi_3 11.5 11.5 13.5 13.5 13.5 13.5 13.5 13.5 16.5 16.5 16.5 16.5 16.5 \
#--zMin 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 --zMax 40 40 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \

plotProfile -m H3K4me3_plots \
--samplesLabel 2C 4C 8C Morula ICM \
--startLabel 'PS' \
--endLabel 'PE' \
-out H3K4me3_plots_profile.pdf

plotHeatmap -m H3K427me3_plots  \
-out H3K27me3_plots.pdf \
--interpolationMethod gaussian \
--whatToShow 'heatmap and colorbar' \
--colorMap Oranges Oranges Oranges Oranges \
--samplesLabel 2C 4C 8C Morula ICM \
--startLabel 'PS' \
--endLabel 'PE' 

plotProfile -m H3K27me3_plots \
--samplesLabel 2C 4C 8C Morula ICM \
--startLabel 'PS' \
--endLabel 'PE' \
-out H3K27me3_plots_profile.pdf

#gunzip *.gz
