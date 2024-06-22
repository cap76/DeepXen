# Multiple Inputs
from keras.utils import plot_model
from keras.models import Model
from keras.layers import Input
from keras.layers import Dense
from keras.layers import Activation
from keras.layers import Flatten
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
Output2 = genfromtxt('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All_GB//allpeaks_labelled_cum_final.s.bed',delimiter='\t',dtype=None)

#Load in the list of genes
onlist = genfromtxt('ChIP_SomiteEndo_All_GB/OnMemFDRtGenes.txt',delimiter='\t',dtype=None)
downlist = genfromtxt('ChIP_SomiteEndo_All_GB/ReprogrammedDowntGenes.txt',delimiter='\t',dtype=None)
offlist = genfromtxt('ChIP_SomiteEndo_All_GB/OffMemFDRtGenes.txt',delimiter='\t',dtype=None)
uplist = genfromtxt('ChIP_SomiteEndo_All_GB/ReprogrammedUptGenes.txt',delimiter='\t',dtype=None)
complist = genfromtxt('ChIP_SomiteEndo_All_GB/ComplementarySettGenes.txt',delimiter='\t',dtype=None)

foxlist = genfromtxt('BW_Endo_All2/FOXA1ChIP.txt',delimiter='\t',dtype=None)


#Generate the Y matrix
Yalt = np.zeros((np.shape(Output2)[0],5))
Xexpa = np.zeros((np.shape(Output2)[0],1,2))

for i in range(0, np.shape(Yalt)[0]):
      Yalt[i,0] = (np.intersect1d(Output2['f3'][i],onlist)).size
      Yalt[i,1] = (np.intersect1d(Output2['f3'][i],downlist)).size
      Yalt[i,2] = (np.intersect1d(Output2['f3'][i],offlist)).size
      Yalt[i,3] = (np.intersect1d(Output2['f3'][i],uplist)).size
      Yalt[i,4] = 1-max(Yalt[i,0:4]) #(np.intersect1d(Output2['f3'][i],complist)).size
      #Xexpa[i,0,2] = (np.intersect1d(Output2['f3'][i],foxlist)).size


#Load in the expression data (do NOT use NT)
#Xexpa = np.zeros((np.shape(Output2)[0],1,4))
#Xexpa[0:np.shape(Output2)[0],0,0] = np.log2(Output2['f4']+1)
#Xexpa[0:np.shape(Output2)[0],0,1] = np.log2(Output2['f5']+1)
#Xexpa[0:np.shape(Output2)[0],0,3] = np.log2(Output2['f6']+1)

#explist = genfromtxt('ChIP_SomiteEcto_All/edger_de_pairing_BIGtable_ectoIVF_somiteDonor_ectoNT_plus_DE.p.csv',delimiter='\t',dtype=None)
#explist = genfromtxt('ChIP_SomiteEndo_All_GB//edger_de_pairing_BIGtable_endoIVF_ectoDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)
explist = genfromtxt('ChIP_SomiteEndo_All_GB/edger_de_pairing_BIGtable_endoIVF_somiteDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)
#Now load the expression for the target tissue: ectoderm/mesoderm?
#for i in range(0, np.shape(Yalt)[0]):
#   arr_ind = np.where(Output2['f3'][i] == explist['f0'])
#   Xexpa[i,0,0] = np.log2(explist['f1'][arr_ind]+1)
#   Xexpa[i,0,1] = np.log2(explist['f2'][arr_ind]+1)

for i in range(0, np.shape(Yalt)[0]):
   arr_ind = np.where(Output2['f3'][i] == explist['f0'])
   Xexpa[i,0,0] = np.log2(explist['f1'][arr_ind].astype('f')+1)
   Xexpa[i,0,1] = np.log2(explist['f2'][arr_ind].astype('f')+1)



np.nan_to_num(Xexpa,copy=False)

#fil=glob.glob("/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All_GB/*.tab")
fil=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All_GB/*.tab'))

#Load in all the expression data
encoder = LabelEncoder()
meh = encoder.fit(Output2['f0'])
Xchr = np.zeros((np.shape(Output2)[0],300,np.shape(fil)[0]))
for k in range(0, np.shape(fil)[0]):
        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
        Xchr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output

#Split into train/val/test sets
trainset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-18]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-17]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-16]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-15]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-14]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-13]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-12]))
valset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-11]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-10]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-9]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-8]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-7]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-6]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-5]))
testset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]))

#Normalise the expression data
Xexp = copy.deepcopy(Xexpa)
Xexp[:,0,0] = ( Xexp[:,0,0] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0])
Xexp[:,0,1] = ( Xexp[:,0,1] - np.nanmean(Xexp[trainset,0,1]) ) / np.nanstd(Xexp[trainset,0,1])
#Xexp[:,0,3] = ( Xexp[:,0,3] - np.nanmean(Xexp[trainset,0,3]) ) / np.nanstd(Xexp[trainset,0,3])

Xchra = copy.deepcopy(Xchr)

#Normalise histone modification data
for k in range(0, np.shape(fil)[0]):
      Xchr[:,0:300,k] = ( Xchr[:,0:300,k] - np.nanmean(Xchr[trainset,0:300,k]) ) / np.nanstd(Xchr[trainset,0:300,k])


from keras.models import load_model
#model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAllGB/Model.h5')
#from keras.models import load_model

#model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/Model.h5')
model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAllGB/Model_seed=2.h5')

#ResultsSomiteEndoAllGB/Model_seed=2.h5

model1 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAllGB/Model_TL_seed=2.h5')


import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
plt.switch_backend('agg')

from sklearn.metrics import roc_curve, precision_recall_curve,auc

fpr_keras, tpr_keras, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())

fpr_keras1, tpr_keras1, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model1.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_keras1, re_keras1, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model1.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())

np.save('ResultsSomiteEndoAllGB/fpr_kears_TL_seed=2.npy', fpr_keras)
np.save('ResultsSomiteEndoAllGB/fpr_kears1_TL_seed=2.npy', fpr_keras1)
np.save('ResultsSomiteEndoAllGB/tpr_kears_TL_seed=2.npy', tpr_keras)
np.save('ResultsSomiteEndoAllGB/tpr_kears1_TL_seed=2.npy', tpr_keras1)
np.save('ResultsSomiteEndoAllGB/pr_kears_TL_seed=2.npy', pr_keras)
np.save('ResultsSomiteEndoAllGB/pr_kears1_TL_seed=2.npy', pr_keras1)
np.save('ResultsSomiteEndoAllGB/re_kears_TL_seed=2.npy', re_keras)
np.save('ResultsSomiteEndoAllGB/re_kears1_TL_seed=2.npy', re_keras1)

classpred = model1.predict([Xchr,Xexp], batch_size=1000)

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

M[0,0]= sum(OnP) #,val=True)
M[0,1]= sum(OnT) #,val=True) #.shape[0]
M[0,2]= sum(OnTP) #,val=True) #.shape[0]
M[0,3]= sum(OnFP) #,val=True) #.shape[0]
M[0,4]= sum(OnTN) #,val=True) #.shape[0]
M[0,5]= sum(OnFN) #,val=True) #.shape[0]
M[0,6]= sum(OffP) #,val=True) #.shape[0]
M[0,7]= sum(OffT) #,val=True) #.shape[0]
M[0,8]= sum(OffTP) #,val=True) #.shape[0]
M[0,9]= sum(OffFP) #,val=True) #.shape[0]
M[0,10]= sum(OnTN) #,val=True) #.shape[0]
M[0,11]= sum(OffFN) #,val=True) #.shape[0]
M[0,12]= sum(RUP) #,val=True) #.shape[0]
M[0,13]= sum(RUT) #,val=True) #.shape[0]
M[0,14]= sum(RUTP) #,val=True) #.shape[0]
M[0,15]= sum(RUFP) #,val=True) #.shape[0]
M[0,16]= sum(RUTN) #,val=True) #.shape[0]
M[0,17]= sum(RUFN) #,val=True) #.shape[0]
M[0,18]= sum(RDP) #,val=True) #.shape[0]
M[0,19]= sum(RDT) #,val=True) #.shape[0]
M[0,20]= sum(RDTP) #,val=True) #.shape[0]
M[0,21]= sum(RDFP) #,val=True) #.shape[0]
M[0,22]= sum(RDTN) #,val=True) #.shape[0]
M[0,23]= sum(RDFN) #,val=True) #.shape[0]
M[0,24]= sum(OtherP) #,val=True) #.shape[0]
M[0,25]= sum(OtherT) #,val=True) #shape[0]
pd.DataFrame(M).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAllGB/Count_TLs_seed=2.csv')
