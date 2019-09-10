%In this folder we have the following files.

Off_RU_Expression.pdf - a summary of the expression profiles of memory genes. Here we indicate Off memory and reprogrammed up expression. 
			log2FC NT-IVF vs Donor (RPKM) as defined in the Xen paper


On_RD_Expression.pdf - a summary of the expression profiles of memory genes. Here we indicate On memory and reprogrammed down expression.
                        log2FC NT-IVF vs Donor (RPKM) as defined in the Xen paper

AUC - ROC curve (evaluated on 1/3 held out data)

PR - Precision-recall curve (evaluated on held out data)

Folder RawActivationsExamples
	In this folder we plot examples of the activation function across around the TSS of various genes for true positives. We do
	this for On, Off, RU and RD. These are for individual genes (or augmented genes i.e., where we duplicate a gene and
	shuffle its position slightly)

ActivationClusters_Memory_minus_Reprogrammed
	In this folder we represent the average activation function AF for all true positives e.g., OFFTP_ac.pdf represents the AF
	for Off memory TPs. We also cluster functions e.g., OFFTP_ac_cl1.pdf represents cluster 1.
	IMPORTANT: these AFs are slightly different to what we've been showing before. Here I indicate the AF relative to the appropriate
	control e.g., for the On memory I indicate AF_On - AF_RD. That is the activation in on memory relatice to RD, so we're seeing
	essentially what the On memory is seeing over and above the RD.
	For Off we compare to RU. For RU we compare to Off. For RD we compare to On.



