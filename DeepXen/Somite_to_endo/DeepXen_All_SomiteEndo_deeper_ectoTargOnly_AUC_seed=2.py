
# Multiple Inputs
from keras.utils import plot_model
from keras.models import Model
from keras.layers import Input
from keras.layers import Dense
from keras.layers import Flatten
from keras.layers import Activation
from keras.layers import Dropout
from keras import optimizers
from keras import regularizers
from keras.layers.convolutional import Conv2D
from keras.layers.convolutional import Conv1D
from keras.layers.pooling import MaxPooling2D
from keras.layers.pooling import MaxPooling1D
from keras.layers.merge import concatenate
from keras.utils import np_utils
from sklearn.preprocessing import LabelEncoder
from numpy import genfromtxt
import numpy as np
import pandas as pd
import glob
import os
import sys
import copy
import matplotlib.pyplot as plt 
from matplotlib import colors as mcolors
plt.switch_backend('agg')

from keras.callbacks import Callback
class GetBest(Callback):
    """Get the best model at the end of training.
	# Arguments
        monitor: quantity to monitor.
        verbose: verbosity mode, 0 or 1.
        mode: one of {auto, min, max}.
            The decision
            to overwrite the current stored weights is made
            based on either the maximization or the
            minimization of the monitored quantity. For `val_acc`,
            this should be `max`, for `val_loss` this should
            be `min`, etc. In `auto` mode, the direction is
            automatically inferred from the name of the monitored quantity.
        period: Interval (number of epochs) between checkpoints.
	# Example
		callbacks = [GetBest(monitor='val_acc', verbose=1, mode='max')]
		mode.fit(X, y, validation_data=(X_eval, Y_eval),
                 callbacks=callbacks)
    """

    def __init__(self, monitor='val_loss', verbose=0,
                 mode='auto', period=1):
        super(GetBest, self).__init__()
        self.monitor = monitor
        self.verbose = verbose
        self.period = period
        self.best_epochs = 0
        self.epochs_since_last_save = 0

        if mode not in ['auto', 'min', 'max']:
            warnings.warn('GetBest mode %s is unknown, '
                          'fallback to auto mode.' % (mode),
                          RuntimeWarning)
            mode = 'auto'

        if mode == 'min':
            self.monitor_op = np.less
            self.best = np.Inf
        elif mode == 'max':
            self.monitor_op = np.greater
            self.best = -np.Inf
        else:
            if 'acc' in self.monitor or self.monitor.startswith('fmeasure'):
                self.monitor_op = np.greater
                self.best = -np.Inf
            else:
                self.monitor_op = np.less
                self.best = np.Inf
                
    def on_train_begin(self, logs=None):
        self.best_weights = self.model.get_weights()

    def on_epoch_end(self, epoch, logs=None):
        logs = logs or {}
        self.epochs_since_last_save += 1
        if self.epochs_since_last_save >= self.period:
            self.epochs_since_last_save = 0
            #filepath = self.filepath.format(epoch=epoch + 1, **logs)
            current = logs.get(self.monitor)
            if current is None:
                warnings.warn('Can pick best model only with %s available, '
                              'skipping.' % (self.monitor), RuntimeWarning)
            else:
                if self.monitor_op(current, self.best):
                    if self.verbose > 0:
                        print('\nEpoch %05d: %s improved from %0.5f to %0.5f,'
                              ' storing weights.'
                              % (epoch + 1, self.monitor, self.best,
                                 current))
                    self.best = current
                    self.best_epochs = epoch + 1
                    self.best_weights = self.model.get_weights()
                else:
                    if self.verbose > 0:
                        print('\nEpoch %05d: %s did not improve' %
                              (epoch + 1, self.monitor))            
                    
    def on_train_end(self, logs=None):
        if self.verbose > 0:
            print('Using epoch %05d with %s: %0.5f' % (self.best_epochs, self.monitor,
                                                       self.best))
        self.model.set_weights(self.best_weights)

np.random.seed(0)

#Process the data. Include this so everyhting is repeatable.
#os.chdir('./ChIP_Endo')
#os.system('./processpeaks_2.sh')
#os.chdir('../BW_Endo')
#os.system('./getextendedpeaks.sh')
#os.chdir('../')


#Load in the raw data (labels are encoded in column 4) and one-hot encode.
#Output2 = genfromtxt('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/allpeaks_labelled_cum_final.s.bed',delimiter='\t',dtype=None)


Output2 = genfromtxt('BW_SomiteEndo_All/allpeaks_labelled_cum_final.se.bed',delimiter='\t',dtype=None)
#EndoNT_On_all.inde.bed',delimiter='\t',dtype=None)

#Load in the list of genes
onlist = genfromtxt('ChIP_SomiteEndo_All/OnMemFDRtGenes.txt',delimiter='\t',dtype=None)
downlist = genfromtxt('ChIP_SomiteEndo_All/ReprogrammedDowntGenes.txt',delimiter='\t',dtype=None)
offlist = genfromtxt('ChIP_SomiteEndo_All/OffMemFDRtGenes.txt',delimiter='\t',dtype=None)
uplist = genfromtxt('ChIP_SomiteEndo_All/ReprogrammedUptGenes.txt',delimiter='\t',dtype=None)
complist = genfromtxt('ChIP_SomiteEndo_All/ComplementarySettGenes.txt',delimiter='\t',dtype=None)

#foxlist = genfromtxt('BW_Endo_All2/FOXA1ChIP.txt',delimiter='\t',dtype=None)

#enhlist = genfromtxt('BW_Endo_All2/enhancers.bed',delimiter='\t',dtype=None)


#Generate the Y matrix
Yalt = np.zeros((np.shape(Output2)[0],5))
Xexpa = np.zeros((np.shape(Output2)[0],1,1))

for i in range(0, np.shape(Yalt)[0]):
      Yalt[i,0] = (np.intersect1d(Output2['f3'][i],onlist)).size
      Yalt[i,1] = (np.intersect1d(Output2['f3'][i],downlist)).size
      Yalt[i,2] = (np.intersect1d(Output2['f3'][i],offlist)).size
      Yalt[i,3] = (np.intersect1d(Output2['f3'][i],uplist)).size
      Yalt[i,4] = 1-max(Yalt[i,0:4]) 

explist = genfromtxt('ChIP_SomiteEndo_All/edger_de_pairing_BIGtable_endoIVF_somiteDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)

#explist = genfromtxt('ChIP_SomiteEndo_All/edger_de_pairing_BIGtable_ectoIVF_somiteDonor_ectoNT_plus_DE.p.csv',delimiter='\t',dtype=None)
#explist = genfromtxt('ChIP_SomiteEcto_All/edger_de_pairing_BIGtable_endoIVF_ectoDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)

for i in range(0, np.shape(Yalt)[0]):
   arr_ind = np.where(Output2['f3'][i] == explist['f0'])
   #Xexpa[i,0,0] = np.log2(explist['f1'][arr_ind].astype('f')+1)
   Xexpa[i,0,0] = np.log2(explist['f2'][arr_ind].astype('f')+1)


np.nan_to_num(Xexpa,copy=False)

#fil=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All/*bwe.tab'))
#DeepXen_All_SomiteEndo_deeper.py
#fil=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All/*bwe.tab'))
#fil1=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All/Endo*bwe.tab')) 
fil2=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All/Ecto*bwe.tab'))
fil3=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All/GC.bwe.tab'))
fil4=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All/*Methylation.bwe.tab'))

fil=fil2+fil3+fil4

#Load in all the expression data
encoder = LabelEncoder()
meh = encoder.fit(Output2['f0'])
Xchr = np.zeros((np.shape(Output2)[0],600,np.shape(fil)[0]))
for k in range(0, np.shape(fil)[0]):
        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
        Xchr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output


#fil2=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/*top*.tab'))

#Outputlab = genfromtxt('BW_Endo_All2/Motif_ap1.bed.bot.bw_oute.bed',delimiter='\t',dtype=None)
#XTF1 = np.zeros((np.shape(Output2)[0],600,np.shape(fil2)[0]))
#XTF2 = np.zeros((np.shape(Output2)[0],600,np.shape(fil2)[0]))


#Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)

#Chromatin mods
#for k in range(0, np.shape(fil)[0]):
#        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
#        Xchr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output

#Better to frame this as an outer join?
#indexss = np.zeros((np.shape(Outputlab)[0],),dtype=np.int8)
#for k in range(0, np.shape(Outputlab)[0]):
#       indexss[k] = np.where((Outputlab['f0'][k]==Output2['f0']) & (Outputlab['f1'][k]==Output2['f1']) & (Outputlab['f2'][k]==Output2['f2']))[0][0]

#Output = genfromtxt('BW_Endo_All2/Motif_ap1.bed.bot.bwe.tab',delimiter='\t',dtype=None, skip_header=3)

#TF binding
#for k in range(0, np.shape(fil2)[0]):
#      Output = genfromtxt(fils2[k],delimiter='\t',dtype=None, skip_header=3)
#      XTF1[indexss.astype(np.int64),:,0] = Output
#      Output = genfromtxt(fils2[k]replace("top","bot"),delimiter='\t',dtype=None, skip_header=3)
#      XTF2[indexss.astype(np.int64),:,0] = Output

#Now merge the three datasets?

#Split into train/val/test sets
trainset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-18]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-17]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-16]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-15]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-14]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-13]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-12]))
valset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-11]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-10]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-9]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-8]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-7]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-6]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-5]))
testset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]))

#Normalise the expression data
Xexp = copy.deepcopy(Xexpa)
Xexp[:,0,0] = ( Xexp[:,0,0] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0])
#Xexp[:,0,1] = ( Xexp[:,0,1] - np.nanmean(Xexp[trainset,0,1]) ) / np.nanstd(Xexp[trainset,0,1])
#Xexp[:,0,3] = ( Xexp[:,0,3] - np.nanmean(Xexp[trainset,0,3]) ) / np.nanstd(Xexp[trainset,0,3])

Xchra = copy.deepcopy(Xchr)

#Normalise histone modification data
for k in range(0, np.shape(fil)[0]):
      Xchr[:,0:600,k] = ( Xchr[:,0:600,k] - np.nanmean(Xchr[trainset,0:600,k]) ) / np.nanstd(Xchr[trainset,0:600,k])

from keras.models import load_model
model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAll/Model_deeper_ectoTargOnly_seed=2.h5')



#from keras.models import load_model
#model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_endoTargOnly.h5')
#plt.switch_backend('agg')
#Scores = np.zeros((1,1))
#Scores[0,0] = model.evaluate([Xchr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
#pd.DataFrame(Scores, columns=['Accuracy']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/PredictionScores_endoTargOnly.csv')
#classpred1 = model.predict([Xchr,Xexp], batch_size=1000)
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
plt.switch_backend('agg')
from sklearn.metrics import roc_curve, precision_recall_curve,auc
#fpr_keras, tpr_keras, thresholds_keras = roc_curve(Yalt.ravel(), model.predict([Xchr,Xexp], batch_size=1000).ravel())
#pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt.ravel(), model.predict([Xchr,Xexp], batch_size=1000).ravel())
#fpr_kerast, tpr_kerast, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
#pr_kerast, re_kerast, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
#np.save('ResultsEndoAll/fpr_kearst_endoTargOnly.npy', fpr_keras)
#np.save('ResultsEndoAll/tpr_kearst_endoTargOnly.npy', tpr_keras)
#np.save('ResultsEndoAll/pr_kearst_endoTargOnly.npy', pr_keras)
#np.save('ResultsEndoAll/re_kearst_endoTargOnly.npy', re_keras)
Y0 = Yalt[valset,:]
P0 = model.predict([Xchr[valset,:,:],Xexp[valset,:,:]])
fpr_keras0_1, tpr_keras0_1, thresholds_kerast = roc_curve(Y0[:,[True,False,False,False,False]], P0[:,[True,False,False,False,False]])
fpr_keras0_2, tpr_keras0_2, thresholds_kerast = roc_curve(Y0[:,[False,True,False,False,False]], P0[:,[False,True,False,False,False]])
fpr_keras0_3, tpr_keras0_3, thresholds_kerast = roc_curve(Y0[:,[False,False,True,False,False]], P0[:,[False,False,True,False,False]])
fpr_keras0_4, tpr_keras0_4, thresholds_kerast = roc_curve(Y0[:,[False,False,False,True,False]], P0[:,[False,False,False,True,False]])
fpr_keras0_5, tpr_keras0_5, thresholds_kerast = roc_curve(Y0[:,[False,False,False,False,True]], P0[:,[False,False,False,False,True]])
np.save('ResultsSomiteEndoAll/tpr_kears0_1_smoendoTargOnly_seed=2.npy', tpr_keras0_1)
np.save('ResultsSomiteEndoAll/tpr_kears0_2_somendoTargOnly_seed=2.npy', tpr_keras0_2)
np.save('ResultsSomiteEndoAll/tpr_kears0_3_somendoTargOnly_seed=2.npy', tpr_keras0_3)
np.save('ResultsSomiteEndoAll/tpr_kears0_4_somendoTargOnly_seed=2.npy', tpr_keras0_4)
np.save('ResultsSomiteEndoAll/tpr_kears0_5_somendoTargOnly_seed=2.npy', tpr_keras0_5)
np.save('ResultsSomiteEndoAll/fpr_kears0_1_somendoTargOnly_seed=2.npy', fpr_keras0_1)
np.save('ResultsSomiteEndoAll/fpr_kears0_2_somendoTargOnly_seed=2.npy', fpr_keras0_2)
np.save('ResultsSomiteEndoAll/fpr_kears0_3_somendoTargOnly_seed=2.npy', fpr_keras0_3)
np.save('ResultsSomiteEndoAll/fpr_kears0_4_somendoTargOnly_seed=2.npy', fpr_keras0_4)
np.save('ResultsSomiteEndoAll/fpr_kears0_5_somendoTargOnly_seed=2.npy', fpr_keras0_5)
pr_keras0_1, re_keras0_1, thresholds_kerast = precision_recall_curve(Y0[:,[True,False,False,False,False]], P0[:,[True,False,False,False,False]])
pr_keras0_2, re_keras0_2, thresholds_kerast = precision_recall_curve(Y0[:,[False,True,False,False,False]], P0[:,[False,True,False,False,False]])
pr_keras0_3, re_keras0_3, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,True,False,False]], P0[:,[False,False,True,False,False]])
pr_keras0_4, re_keras0_4, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,False,True,False]], P0[:,[False,False,False,True,False]])
pr_keras0_5, re_keras0_5, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,False,False,True]], P0[:,[False,False,False,False,True]])

np.save('ResultsSomiteEndoAll/pr_kears0_1_somendoTargOnly_seed=2.npy', pr_keras0_1)
np.save('ResultsSomiteEndoAll/pr_kears0_2_somendoTargOnly_seed=2.npy', pr_keras0_2)
np.save('ResultsSomiteEndoAll/pr_kears0_3_somendoTargOnly_seed=2.npy', pr_keras0_3)
np.save('ResultsSomiteEndoAll/pr_kears0_4_somendoTargOnly_seed=2.npy', pr_keras0_4)
np.save('ResultsSomiteEndoAll/pr_kears0_5_somendoTargOnly_seed=2.npy', pr_keras0_5)

np.save('ResultsSomiteEndoAll/re_kears0_1_somendoTargOnly_seed=2.npy', re_keras0_1)
np.save('ResultsSomiteEndoAll/re_kears0_2_somendoTargOnly_seed=2.npy', re_keras0_2)
np.save('ResultsSomiteEndoAll/re_kears0_3_somendoTargOnly_seed=2.npy', re_keras0_3)
np.save('ResultsSomiteEndoAll/re_kears0_4_somendoTargOnly_seed=2.npy', re_keras0_4)
np.save('ResultsSomiteEndoAll/re_kears0_5_somendoTargOnly_seed=2.npy', re_keras0_5)







#model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAll/Model_deeper_ectoOnly.h5')
plt.switch_backend('agg')
Scores = np.zeros((1,1))
Scores[0,0] = model.evaluate([Xchr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
pd.DataFrame(Scores, columns=['Accuracy']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAll/PredictionScores_ectoTargOnly_seed=2.csv')
classpred1 = model.predict([Xchr,Xexp], batch_size=1000)
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
plt.switch_backend('agg')
from sklearn.metrics import roc_curve, precision_recall_curve,auc
fpr_keras, tpr_keras, thresholds_keras = roc_curve(Yalt.ravel(), model.predict([Xchr,Xexp], batch_size=1000).ravel())
pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt.ravel(), model.predict([Xchr,Xexp], batch_size=1000).ravel())
fpr_kerast, tpr_kerast, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_kerast, re_kerast, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
np.save('ResultsSomiteEndoAll/fpr_kearst_ectoTargOnly_seed=2.npy', fpr_keras)
np.save('ResultsSomiteEndoAll/tpr_kearst_ectoTargOnly_seed=2.npy', tpr_keras)
np.save('ResultsSomiteEndoAll/pr_kearst_ectoTargOnly_seed=2.npy', pr_keras)
np.save('ResultsSomiteEndoAll/re_kearst_ectoTargOnly_seed=2.npy', re_keras)
plt.figure()
plt.plot(re_kerast, pr_kerast, label='Combined') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig('ResultsSomiteEndoAll/PR_ectoTargOnly_seed=2.pdf')
plt.draw()
plt.figure()
plt.plot(fpr_kerast, tpr_kerast, label='Combined') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot([0, 1], [0, 1], 'k--', label='Random')
plt.plot([0,0, 1], [0, 1, 1], 'k-')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig('ResultsSomiteEndoAll/AUC_ectoTargOnly_seed=2.pdf')
plt.draw()
classpred = classpred1
OnP = ((classpred[:,0]>np.max(classpred[:,1:5],axis=1)))
OnT = (Yalt[:,0]==1)
OnTP = ((classpred[:,0]>np.max(classpred[:,1:5],axis=1)) & (Yalt[:,0]==1))
OnFP = ((classpred[:,0]>np.max(classpred[:,1:5],axis=1)) & (Yalt[:,0]==0))
OnTN = ((classpred[:,0]<np.max(classpred[:,1:5],axis=1)) & (Yalt[:,0]==0))
OnFN = ((classpred[:,0]<np.max(classpred[:,1:5],axis=1)) & (Yalt[:,0]==1))
OffT = (Yalt[:,2]==1)
OffP = ((classpred[:,2]>np.max(classpred[:,[True,True,False,True,True]],axis=1) ))
OffTP = ((classpred[:,2]>np.max(classpred[:,[True,True,False,True,True]],axis=1)) & (Yalt[:,2]==1))
OffFP = ((classpred[:,2]>np.max(classpred[:,[True,True,False,True,True]],axis=1)) & (Yalt[:,2]==0))
OffTN = ((classpred[:,2]<np.max(classpred[:,[True,True,False,True,True]],axis=1)) & (Yalt[:,2]==0))
OffFN = ((classpred[:,2]<np.max(classpred[:,[True,True,False,True,True]],axis=1)) & (Yalt[:,2]==1))
RUT = (Yalt[:,3]==1)
RUP = ((classpred[:,3]>np.max(classpred[:,[True,True,True,False,True]],axis=1)))
RUTP = ((classpred[:,3]>np.max(classpred[:,[True,True,True,False,True]],axis=1)) & (Yalt[:,3]==1))
RUFP = ((classpred[:,3]>np.max(classpred[:,[True,True,True,False,True]],axis=1)) & (Yalt[:,3]==0))
RUTN = ((classpred[:,3]<np.max(classpred[:,[True,True,True,False,True]],axis=1)) & (Yalt[:,3]==0))
RUFN = ((classpred[:,3]<np.max(classpred[:,[True,True,True,False,True]],axis=1)) & (Yalt[:,3]==1))
RDT = (Yalt[:,1]==1)
RDP = ((classpred[:,1]>np.max(classpred[:,[True,False,True,True,True]],axis=1)))
RDTP = ((classpred[:,1]>np.max(classpred[:,[True,False,True,True,True]],axis=1)) & (Yalt[:,1]==1))
RDFP = ((classpred[:,1]>np.max(classpred[:,[True,False,True,True,True]],axis=1)) & (Yalt[:,1]==0))
RDTN = ((classpred[:,1]<np.max(classpred[:,[True,False,True,True,True]],axis=1)) & (Yalt[:,1]==0))
RDFN = ((classpred[:,1]<np.max(classpred[:,[True,False,True,True,True]],axis=1)) & (Yalt[:,1]==1))
OtherP = ((classpred[:,4]>np.max(classpred[:,[True,True,True,True,False]],axis=1)))
OtherT = (Yalt[:,4]==1)
M = np.zeros((1,26))
M[0,0]= sum(OnP)
M[0,1]= sum(OnT)
M[0,2]= sum(OnTP)
M[0,3]= sum(OnFP)
M[0,4]= sum(OnTN)
M[0,5]= sum(OnFN)
M[0,6]= sum(OffP)
M[0,7]= sum(OffT)
M[0,8]= sum(OffTP)
M[0,9]= sum(OffFP)
M[0,10]= sum(OnTN)
M[0,11]= sum(OffFN)
M[0,12]= sum(RUP)
M[0,13]= sum(RUT)
M[0,14]= sum(RUTP)
M[0,15]= sum(RUFP)
M[0,16]= sum(RUTN)
M[0,17]= sum(RUFN)
M[0,18]= sum(RDP)
M[0,19]= sum(RDT)
M[0,20]= sum(RDTP)
M[0,21]= sum(RDFP)
M[0,22]= sum(RDTN)
M[0,23]= sum(RDFN)
M[0,24]= sum(OtherP)
M[0,25]= sum(OtherT)
pd.DataFrame(M).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAll/Counts_ectoTargOnly_seed=2.csv')

