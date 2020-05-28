
#Analysis for DeepXen as follows:

EndoFull: the full DeepXen model trained on chromatin marks and expression values
EndoChr - model for chromatin marks only
EndoTL - trasfer learning using model from SomEndo to predict Endo data

#Corresponding folders for the SomEndo datasets
EndoFullSomEndo 
EndoChrSomEndo
SomEndoTL


#In each folder will be a number of pdf files:

RUTP_ac.pdf, RUFP_ac.pdf etc are the based activations i.e., what marks are contributing
most to action of TP

RDTP_Deltac.pdf, etc are the delta activations for the corresponding pairs
 i.e., RD activation - On activatio

These are the overal actions, I have also clustered the activations into 5 groups
indicated by clN where N is an arbitrary cluster number.

There is also a second run indicated by _seed=2 for each activation plot. The clusers
will not necessarily correspond. But the overall acivations should be comparable.
Note that there may be some difference in the scales, and we may need to manually set
these for final plots. 

There are also corresponding set of plots for the BG datasets.


#Models for the GB

#SomEndoGB
ResultsSomiteEndoAllGB/Model_seed=2.h5 - full
Model.h5 - full
Model_pos_seed=2.h5
Model_pos.h5
Model_neg_seed=2.h5
Model_neg.h5

#Endo
ResultsEndoAllGB/Model_seed=2.h5
Model_seed=2.h5 - full
Model_contneg_seed=2.h5
Model_contneg.h5
Model_cont_seed=2.h5
Model_cont.h5

#Models for TSS

