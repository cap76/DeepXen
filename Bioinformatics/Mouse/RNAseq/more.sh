#mkdir ./SRR7007738
#cd ./SRR7007738
#STAR --genomeDir /mnt/scratch/surani/cap76/Lena/RNAseq/STARgenome/ --sjdbGTFfile /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Annotation/Genes/genes.gtf --sjdbOverhang 100 --readFilesIn ../SRR7007738_1_val_1.fq ../SRR7007738_2_val_2.fq
#htseq-count -m intersection-nonempty -i gene_name Aligned.out.sam /mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Annotation/Genes/genes.gtf > SRR7007738.counts

#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR700/004/SRR7007704/SRR7007704_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR700/004/SRR7007704/SRR7007704_2.fastq.gz
#gunzip SRR7007704_1.fastq.gz
#gunzip SRR7007704_2.fastq.gz
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse/RNAseq --paired SRR7007704_1.fastq  SRR7007704_2.fastq 


mkdir SRR7007704
cd SRR7007704
STAR --genomeDir /mnt/scratch/surani/cap76/Lena/RNAseq/STARgenome/ --sjdbGTFfile /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Annotation/Genes/genes.gtf --sjdbOverhang 100 --readFilesIn ../SRR7007704_1_val_1.fq ../SRR7007704_2_val_2.fq
htseq-count -m intersection-nonempty -i gene_name Aligned.out.sam /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Annotation/Genes/genes.gtf > SRR7007704.counts

