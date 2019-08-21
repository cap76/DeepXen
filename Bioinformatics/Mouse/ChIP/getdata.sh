#GSM1483900	MEF.H3.MNase-ChIP-Seq
#GSM1483901	MEF.H3K4me3.ChIP-Seq
#GSM1483902	MEF.H3K27me3.ChIP-Seq
#GSM1483903	MEF.H3K9me3.ChIP-Seq
#GSM1483901
#GSM1816298	MEF_rep.H3.MNase-ChIP-Seq
#GSM1816299	MEF_rep.H3K4me3.ChIP-Seq
#GSM1816300	MEF_rep.H3K27me3.ChIP-Seq

#GSM2417085MEF_H3K4me1_ChIP-Seq
#GSM2417093MEF_H3K27ac_ChIP-Seq
#GSM2417121MEF_Inputnative_MNase_ChIP-Seq

#GSM3399945	H3K4me3 ChIP-seq MEF
#GSM3399946	H3K27ac ChIP-seq MEF
#GSM3399947	H3K27me3 ChIP-seq MEF

#GSM2629919	ChIP-seq_H3K4me3_MEF_rep1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/000/SRR5572610/SRR5572610.fastq.gz
#GSM2629920	ChIP-seq_H3K4me3_MEF_rep2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/001/SRR5572611/SRR5572611.fastq.gz
#GSM2629921	ChIP-seq_H3K4me1_MEF_rep1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/002/SRR5572612/SRR5572612.fastq.gz
#GSM2629922	ChIP-seq_H3K4me1_MEF_rep2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/003/SRR5572613/SRR5572613.fastq.gz
#GSM2629923	ChIP-seq_H3K27ac_MEF_rep1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/004/SRR5572614/SRR5572614.fastq.gz
#GSM2629924	ChIP-seq_H3K27ac_MEF_rep2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/005/SRR5572615/SRR5572615.fastq.gz
#GSM2629925	ChIP-seq_H3K27me3_MEF_rep1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/006/SRR5572616/SRR5572616.fastq.gz
#GSM2629926	ChIP-seq_H3K27me3_MEF_rep2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/007/SRR5572617/SRR5572617.fastq.gz
#GSM2629927	ChIP-seq_H3K9me3_MEF_rep1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/008/SRR5572618/SRR5572618.fastq.gz
#GSM2629928	ChIP-seq_H3K9me3_MEF_rep2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/009/SRR5572619/SRR5572619.fastq.gz
#GSM2629933	ChIP-seq_INPUT_MEF_rep1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/004/SRR5572624/SRR5572624.fastq.gz
#GSM2629934	ChIP-seq_INPUT_MEF_rep2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/005/SRR5572625/SRR5572625.fastq.gz

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

#GSM1438699	h2bub_WT_ChIPSeq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518433/SRR1518433.fastq.gz
#GSM1438701	input_WT_ChIPSeq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/005/SRR1518435/SRR1518435.fastq.gz
#GSM1833646	H3K4me3_WT_ChIPSeq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR213/002/SRR2131052/SRR2131052.fastq.gz
#GSM1833648	H3K79me3_WT_ChIPSeq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR213/004/SRR2131054/SRR2131054.fastq.gz


#SRR5572610  SRR5572624
#SRR5572611 SRR5572625
#SRR5572612 SRR5572624
#SRR5572613 SRR5572625
#SRR5572614 SRR5572624
#SRR5572615 SRR5572625
#SRR5572616 SRR5572624
#SRR5572617 SRR5572625
#SRR5572618 SRR5572624
#SRR5572619 SRR5572625
#SRR1518433 SRR1518435
#SRR2131052 SRR1518435
#SRR2131054 SRR1518435





~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572610.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572611.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572612.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572613.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572614.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572615.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572616.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572617.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572618.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572619.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572624.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR5572625.fastq

~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR1518433.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR1518435.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR2131052.fastq
~/./trim_galore --fastqc --three_prime_clip_R1 1 --clip_R1 4 --output_dir /mnt/scratch/gurdon/cap76/ForMing/Chris SRR2131054.fastq

