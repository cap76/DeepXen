
#rm SRR2131054.fastq
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR213/004/SRR2131054/SRR2131054.fastq.gz

#gunzip SRR2131054.fastq.gz
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572610.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572611.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572612.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572613.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572614.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572615.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572616.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572617.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572618.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572619.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572624.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572625.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR1518433.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR1518435.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR2131052.fastq
#~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/DeepXen/Mouse/ChIP SRR2131054.fastq

#cp /mnt/scratch/gurdon/cap76/ForMing/Chris/SRR55726* ./
#SRR2131054_trimmed_unique.bam

#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572610_trimmed.fq -S SRR5572610_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572611_trimmed.fq -S SRR5572611_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572612_trimmed.fq -S SRR5572612_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572613_trimmed.fq -S SRR5572613_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572614_trimmed.fq -S SRR5572614_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572615_trimmed.fq -S SRR5572615_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572616_trimmed.fq -S SRR5572616_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572617_trimmed.fq -S SRR5572617_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572618_trimmed.fq -S SRR5572618_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572619_trimmed.fq -S SRR5572619_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572624_trimmed.fq -S SRR5572624_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR5572625_trimmed.fq -S SRR5572625_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR1518433_trimmed.fq -S SRR1518433_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR1518435_trimmed.fq -S SRR1518435_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR2131052_trimmed.fq -S SRR2131052_trimmed.sam
#bowtie2 --local --very-sensitive-local -x /mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Sequence/Bowtie2Index/genome -U SRR2131054_trimmed.fq -S SRR2131054_trimmed.sam




#samtools view -bS SRR2131054_trimmed.sam | samtools sort - SRR2131054_trimmed_sorted
#samtools index SRR2131054_trimmed_sorted.bam  SRR21310542_trimmed_sorted.bai
#samtools view -b -q 5 SRR2131054_trimmed_sorted.bam > SRR2131054_trimmed_unique.bam
#samtools index SRR2131054_trimmed_unique.bam SRR2131054_trimmed_unique.bai

#samtools view -bS SRR2131052_trimmed.sam | samtools sort - SRR2131052_trimmed_sorted
#samtools index SRR2131052_trimmed_sorted.bam  SRR2131052_trimmed_sorted.bai
#samtools view -b -q 5 SRR2131052_trimmed_sorted.bam > SRR2131052_trimmed_unique.bam
#samtools index SRR2131052_trimmed_unique.bam SRR2131052_trimmed_unique.bai


#samtools view -bS SRR1518435_trimmed.sam | samtools sort - SRR1518435_trimmed_sorted
#samtools index SRR1518435_trimmed_sorted.bam  SRR1518435_trimmed_sorted.bai
#samtools view -b -q 5 SRR1518435_trimmed_sorted.bam > SRR1518435_trimmed_unique.bam
#samtools index SRR1518435_trimmed_unique.bam SRR1518435_trimmed_unique.bai



#samtools view -bS SRR1518433_trimmed.sam | samtools sort - SRR1518433_trimmed_sorted
#samtools index SRR1518433_trimmed_sorted.bam  SRR1518433_trimmed_sorted.bai
#samtools view -b -q 5 SRR1518433_trimmed_sorted.bam > SRR1518433_trimmed_unique.bam
#samtools index SRR1518433_trimmed_unique.bam SRR1518433_trimmed_unique.bai


#samtools view -bS SRR5572625_trimmed.sam | samtools sort - SRR5572625_trimmed_sorted
#samtools index SRR5572625_trimmed_sorted.bam  SRR5572625_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572625_trimmed_sorted.bam > SRR5572625_trimmed_unique.bam
#samtools index SRR5572625_trimmed_unique.bam SRR5572625_trimmed_unique.bai

#samtools view -bS SRR5572624_trimmed.sam | samtools sort - SRR5572624_trimmed_sorted
#samtools index SRR5572624_trimmed_sorted.bam  SRR5572624_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572624_trimmed_sorted.bam > SRR5572624_trimmed_unique.bam
#samtools index SRR5572624_trimmed_unique.bam SRR5572624_trimmed_unique.bai

#samtools view -bS SRR5572610_trimmed.sam | samtools sort - SRR5572610_trimmed_sorted
#samtools index SRR5572610_trimmed_sorted.bam  SRR5572610_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572610_trimmed_sorted.bam > SRR5572610_trimmed_unique.bam
#samtools index SRR5572610_trimmed_unique.bam SRR5572610_trimmed_unique.bai

#samtools view -bS SRR5572611_trimmed.sam | samtools sort - SRR5572611_trimmed_sorted
#samtools index SRR5572611_trimmed_sorted.bam  SRR5572611_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572611_trimmed_sorted.bam > SRR5572611_trimmed_unique.bam
#samtools index SRR5572611_trimmed_unique.bam SRR5572611_trimmed_unique.bai

#samtools view -bS SRR5572612_trimmed.sam | samtools sort - SRR5572612_trimmed_sorted
#samtools index SRR5572612_trimmed_sorted.bam  SRR5572612_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572612_trimmed_sorted.bam > SRR5572612_trimmed_unique.bam
#samtools index SRR5572612_trimmed_unique.bam SRR5572612_trimmed_unique.bai

#samtools view -bS SRR5572613_trimmed.sam | samtools sort - SRR5572613_trimmed_sorted
#samtools index SRR5572613_trimmed_sorted.bam  SRR5572613_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572613_trimmed_sorted.bam > SRR5572613_trimmed_unique.bam
#samtools index SRR5572613_trimmed_unique.bam SRR5572613_trimmed_unique.bai

#samtools view -bS SRR5572614_trimmed.sam | samtools sort - SRR5572614_trimmed_sorted
#samtools index SRR5572614_trimmed_sorted.bam  SRR5572614_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572614_trimmed_sorted.bam > SRR5572614_trimmed_unique.bam
#samtools index SRR5572614_trimmed_unique.bam SRR5572614_trimmed_unique.bai

#samtools view -bS SRR5572615_trimmed.sam | samtools sort - SRR5572615_trimmed_sorted
#samtools index SRR5572615_trimmed_sorted.bam  SRR5572615_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572615_trimmed_sorted.bam > SRR5572615_trimmed_unique.bam
#samtools index SRR5572615_trimmed_unique.bam SRR5572615_trimmed_unique.bai

#samtools view -bS SRR5572616_trimmed.sam | samtools sort - SRR5572616_trimmed_sorted
#samtools index SRR5572616_trimmed_sorted.bam  SRR5572616_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572616_trimmed_sorted.bam > SRR5572616_trimmed_unique.bam
#samtools index SRR5572616_trimmed_unique.bam SRR5572616_trimmed_unique.bai

#samtools view -bS SRR5572617_trimmed.sam | samtools sort - SRR5572617_trimmed_sorted
#samtools index SRR5572617_trimmed_sorted.bam  SRR5572617_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572617_trimmed_sorted.bam > SRR5572617_trimmed_unique.bam
#samtools index SRR5572617_trimmed_unique.bam SRR5572617_trimmed_unique.bai

#samtools view -bS SRR5572618_trimmed.sam | samtools sort - SRR5572618_trimmed_sorted
#samtools index SRR5572618_trimmed_sorted.bam  SRR5572618_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572618_trimmed_sorted.bam > SRR5572618_trimmed_unique.bam
#samtools index SRR5572618_trimmed_unique.bam SRR5572618_trimmed_unique.bai

#samtools view -bS SRR5572619_trimmed.sam | samtools sort - SRR5572619_trimmed_sorted
#samtools index SRR5572619_trimmed_sorted.bam  SRR5572619_trimmed_sorted.bai
#samtools view -b -q 5 SRR5572619_trimmed_sorted.bam > SRR5572619_trimmed_unique.bam
#samtools index SRR5572619_trimmed_unique.bam SRR5572619_trimmed_unique.bai

#samtools view -bS ERR1328331_trimmed.sam | samtools sort - ERR1328331_trimmed_sorted
#samtools index ERR1328331_trimmed_sorted.bam  ERR1328331_trimmed_sorted.bai

#samtools view -bS ERR1328332_trimmed.sam | samtools sort - ERR1328332_trimmed_sorted
#samtools index ERR1328332_trimmed_sorted.bam  ERR1328332_trimmed_sorted.bai

#samtools view -bS ERR1328333_trimmed.sam | samtools sort - ERR1328333_trimmed_sorted
#samtools index ERR1328333_trimmed_sorted.bam  ERR1328333_trimmed_sorted.bai

#samtools view -b -q 5 ERR1328330_trimmed_sorted.bam > ERR1328330_trimmed_sorted_unique.bam
#samtools view -b -q 5 ERR1328331_trimmed_sorted.bam > ERR1328331_trimmed_sorted_unique.bam
#samtools view -b -q 5 ERR1328332_trimmed_sorted.bam > ERR1328332_trimmed_sorted_unique.bam
#samtools view -b -q 5 ERR1328333_trimmed_sorted.bam > ERR1328333_trimmed_sorted_unique.bam

#samtools index ERR1328330_trimmed_sorted_unique.bam ERR1328330_trimmed_sorted_unique.bai
#samtools index ERR1328331_trimmed_sorted_unique.bam ERR1328331_trimmed_sorted_unique.bai
#samtools index ERR1328332_trimmed_sorted_unique.bam ERR1328332_trimmed_sorted_unique.bai
#samtools index ERR1328333_trimmed_sorted_unique.bam ERR1328333_trimmed_sorted_unique.bai

#macs2 callpeak -t ERR1328330_trimmed_sorted_unique.bam -c ERR1328333_trimmed_sorted_unique.bam -g mm -n Rep1 -q 0.01 --keep-dup all --nomodel 
#macs2 callpeak -t ERR1328331_trimmed_sorted_unique.bam -c ERR1328333_trimmed_sorted_unique.bam -g mm -n Rep2 -q 0.01 --keep-dup all --nomodel 
#macs2 callpeak -t ERR1328332_trimmed_sorted_unique.bam -c ERR1328333_trimmed_sorted_unique.bam -g mm -n Rep3 -q 0.01 --keep-dup all --nomodel 

#bamCompare -b1  ERR1328330_trimmed_sorted_unique.bam -b2 ERR1328333_trimmed_sorted_unique.bam --binSize 10 --extendReads 150 --normalizeUsing CPM  --operation subtract -o Rep1.bw
#bamCompare -b1  ERR1328331_trimmed_sorted_unique.bam -b2 ERR1328333_trimmed_sorted_unique.bam --binSize 10 --extendReads 150 --normalizeUsing CPM  --operation subtract -o Rep2.bw

#bamCompare -b1  ERR1328332_trimmed_sorted_unique.bam -b2 ERR1328333_trimmed_sorted_unique.bam --binSize 10 --extendReads 150 --normalizeUsing CPM  --operation subtract -o Rep3.bw


#SRR5572610  SRR5572624 #H3K4me3
#SRR5572611 SRR5572625 #
#SRR5572612 SRR5572624 #H3K4me1
#SRR5572613 SRR5572625 
#SRR5572614 SRR5572624 #H3K27ac
#SRR5572615 SRR5572625
#SRR5572616 SRR5572624 #H3K27me3
#SRR5572617 SRR5572625
#SRR5572618 SRR5572624 #H3K9me3
#SRR5572619 SRR5572625

#SRR1518433 SRR1518435 #H2Aub
#SRR2131052 SRR1518435 #H3K4me3
#SRR2131054 SRR1518435 #H3K79me3

#samtools merge H3K4me3.bam SRR5572610_trimmed_unique.bam SRR5572611_trimmed_unique.bam
#samtools merge H3K4me1.bam SRR5572612_trimmed_unique.bam SRR5572613_trimmed_unique.bam
#samtools merge H3K27ac.bam SRR5572614_trimmed_unique.bam SRR5572615_trimmed_unique.bam
#samtools merge H3K27me3.bam SRR5572616_trimmed_unique.bam SRR5572617_trimmed_unique.bam
#samtools merge H3K9me3.bam SRR5572618_trimmed_unique.bam SRR5572619_trimmed_unique.bam
#samtools merge Input.bam SRR5572624_trimmed_unique.bam SRR5572625_trimmed_unique.bam

#samtools index H3K4me3.bam H3K4me3.bai
#samtools index H3K4me1.bam H3K4me1.bai
#samtools index H3K27ac.bam H3K27ac.bai
#samtools index H3K27me3.bam H3K27me3.bai
#samtools index H3K9me3.bam H3K9me3.bai
#samtools index Input.bam Input.bai#

#samtools index SRR2131054_trimmed_unique.bam SRR2131054_trimmed_unique.bai
#bamCompare -b1 H3K4me3.bam -b2 Input.bam --binSize 10 --extendReads 150 --normalizeUsing CPM  --operation subtract -o H3K4me3.bw
#bamCompare -b1 H3K4me1.bam -b2 Input.bam --binSize 10 --extendReads 150 --normalizeUsing CPM  --operation subtract -o H3K4me1.bw
#bamCompare -b1 H3K27ac.bam -b2 Input.bam --binSize 10 --extendReads 150 --normalizeUsing CPM  --operation subtract -o H3K27ac.bw
#bamCompare -b1 H3K27me3.bam -b2 Input.bam --binSize 10 --extendReads 150 --normalizeUsing CPM  --operation subtract -o H3K27me3.bw
#bamCompare -b1 H3K9me3.bam -b2 Input.bam --binSize 10 --extendReads 150 --normalizeUsing CPM  --operation subtract -o H3K9me3.bw

#bamCompare -b1 SRR1518433_trimmed_unique.bam -b2 SRR1518435_trimmed_unique.bam --binSize 10 --extendReads 150 --normalizeUsing CPM  --operation subtract -o H2Aub.bw
#bamCompare -b1 SRR2131052_trimmed_unique.bam -b2 SRR1518435_trimmed_unique.bam --binSize 10 --extendReads 150 --normalizeUsing CPM  --operation subtract -o H2Aub.bw
bamCompare -b1 SRR2131054_trimmed_unique.bam -b2 SRR1518435_trimmed_unique.bam --binSize 10 --extendReads 150 --normalizeUsing CPM  --operation subtract -o H3K79me3.bw



#bamCompare -b1 SRR5572610_trimmed_unique.bam  SRR5572624
#bamCompare -b1 SRR5572611_trimmed_unique.bam SRR5572625
#bamCompare -b1 SRR5572612_trimmed_unique.bam SRR5572624
#bamCompare -b1 SRR5572613_trimmed_unique.bam SRR5572625
#bamCompare -b1 SRR5572614_trimmed_unique.bam SRR5572624
#bamCompare -b1 SRR5572615_trimmed_unique.bam SRR5572625
#bamCompare -b1 SRR5572616 SRR5572624
#bamCompare -b1 SRR5572617 SRR5572625
#bamCompare -b1 SRR5572618 SRR5572624
#bamCompare -b1 SRR5572619 SRR5572625
#bamCompare -b1 SRR1518433 SRR1518435
#bamCompare -b1 SRR2131052 SRR1518435
#bamCompare -b1 SRR2131054 SRR1518435


#GSM1483900     MEF.H3.MNase-ChIP-Seq
#GSM1483901     MEF.H3K4me3.ChIP-Seq
#GSM1483902     MEF.H3K27me3.ChIP-Seq
#GSM1483903     MEF.H3K9me3.ChIP-SeqGSM1483901

#GSM1816298     MEF_rep.H3.MNase-ChIP-Seq
#GSM1816299     MEF_rep.H3K4me3.ChIP-Seq
#GSM1816300     MEF_rep.H3K27me3.ChIP-Seq


#GSM2417085MEF_H3K4me1_ChIP-Seq
#GSM2417093MEF_H3K27ac_ChIP-Seq
#GSM2417121MEF_Inputnative_MNase_ChIP-Seq

#GSM3399945     H3K4me3 ChIP-seq MEF
#GSM3399946     H3K27ac ChIP-seq MEF
#GSM3399947     H3K27me3 ChIP-seq MEF
