~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse/RNAseq --paired SRR7007729_1.fastq  SRR7007729_2.fastq 
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse/RNAseq --paired SRR7007730_1.fastq  SRR7007730_2.fastq 

mkdir ./SRR7007729
cd ./SRR7007729
STAR --genomeDir /mnt/scratch/surani/cap76/Lena/RNAseq/STARgenome/ --sjdbGTFfile /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Annotation/Genes/genes.gtf --sjdbOverhang 100 --readFilesIn ../SRR7007729_1_val_1.fq ../SRR7007729_2_val_2.fq
htseq-count -m intersection-nonempty -i gene_name Aligned.out.sam /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Annotation/Genes/genes.gtf > SRR7007729.counts

mkdir ../SRR7007730
cd ../SRR7007730
STAR --genomeDir /mnt/scratch/surani/cap76/Lena/RNAseq/STARgenome/ --sjdbGTFfile /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Annotation/Genes/genes.gtf --sjdbOverhang 100 --readFilesIn ../SRR7007730_1_val_1.fq ../SRR7007730_2_val_2.fq
htseq-count -m intersection-nonempty -i gene_name Aligned.out.sam /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Annotation/Genes/genes.gtf > SRR7007730.counts

