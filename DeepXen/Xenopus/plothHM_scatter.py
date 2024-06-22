#Run feature extraction
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
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
from deepexplain.tensorflow import DeepExplain
from keras import backend as K
from sklearn.metrics import roc_curve, precision_recall_curve,auc
from sklearn.cluster import KMeans

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
Output2 = genfromtxt('BW_Endo_All2/allpeaks_labelled_cum_final.se.bed',delimiter='\t',dtype=None)
#EndoNT_On_all.inde.bed',delimiter='\t',dtype=None)

#Load in the list of genes
onlist = genfromtxt('ChIP_Endo_All/OnMemFDRtGenes.txt',delimiter='\t',dtype=None)
downlist = genfromtxt('ChIP_Endo_All/ReprogrammedDowntGenes.txt',delimiter='\t',dtype=None)
offlist = genfromtxt('ChIP_Endo_All/OffMemFDRtGenes.txt',delimiter='\t',dtype=None)
uplist = genfromtxt('ChIP_Endo_All/ReprogrammedUptGenes.txt',delimiter='\t',dtype=None)
complist = genfromtxt('ChIP_Endo_All/ComplementarySettGenes.txt',delimiter='\t',dtype=None)

#Generate the Y matrix
Yalt = np.zeros((np.shape(Output2)[0],5))
Xexpa = np.zeros((np.shape(Output2)[0],1,2))
Xrawexp = np.zeros((np.shape(Output2)[0],7))

for i in range(0, np.shape(Yalt)[0]):
      Yalt[i,0] = (np.intersect1d(Output2['f3'][i],onlist)).size
      Yalt[i,1] = (np.intersect1d(Output2['f3'][i],downlist)).size
      Yalt[i,2] = (np.intersect1d(Output2['f3'][i],offlist)).size
      Yalt[i,3] = (np.intersect1d(Output2['f3'][i],uplist)).size
      Yalt[i,4] = 1-max(Yalt[i,0:4]) #(np.intersect1d(Output2['f3'][i],complist)).size
      #Xexpa[i,0,2] = (np.intersect1d(Output2['f3'][i],foxlist)).size

#explist = genfromtxt('edger_de_pairing_BIGtable_endoIVF_ectoDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)
#explist1 = genfromtxt('edger_de_pairing_BIGtable_endoIVF_ectoDonor_endoNT_plus_DE.rmnan.csv',delimiter='\t',dtype=None)

explist = genfromtxt('ChIP_Endo_All/edger_de_pairing_BIGtable_endoIVF_ectoDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)
explist1 = genfromtxt('ChIP_Endo_All/edger_de_pairing_BIGtable_endoIVF_ectoDonor_endoNT_plus_DE.rmnan.csv',delimiter='\t',dtype=None)

#fil=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/*bwe.tab'))



#Now load the expression for the target tissue: ectoderm/mesoderm?
for i in range(0, np.shape(Yalt)[0]):
   arr_ind = np.where(Output2['f3'][i] == explist['f0'])
   Xexpa[i,0,0] = np.log2(explist['f1'][arr_ind].astype('f')+1)
   Xexpa[i,0,1] = np.log2(explist['f2'][arr_ind].astype('f')+1)
   Xrawexp[i,0] = (explist1['f1'][arr_ind].astype('f'))
   Xrawexp[i,1] = (explist1['f2'][arr_ind].astype('f'))
   Xrawexp[i,2] = (explist1['f3'][arr_ind].astype('f'))
   Xrawexp[i,3] = (explist1['f4'][arr_ind].astype('f'))
   Xrawexp[i,4] = (explist1['f5'][arr_ind].astype('f'))
   Xrawexp[i,5] = (explist1['f6'][arr_ind].astype('f'))
   Xrawexp[i,6] = (explist1['f7'][arr_ind].astype('f'))

np.nan_to_num(Xexpa,copy=False)

fil=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/*bwe.tab'))

#fil=sorted(glob.glob('/Users/christopherpenfold/Desktop/Code/deepexplain/TSS_deeper/Endoderm/*.tab'))

#Load in all the expression data
encoder = LabelEncoder()
meh = encoder.fit(Output2['f0'])
Xchr = np.zeros((np.shape(Output2)[0],600,np.shape(fil)[0]))
for k in range(0, np.shape(fil)[0]):
        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
        Xchr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output

#Split into train/val/test sets
trainset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-18]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-17]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-16]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-15]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-14]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-13]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-12]))
valset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-11]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-10]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-9]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-8]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-7]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-6]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-5]))
testset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]))


#testset1=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]))
#testset2=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]))
#testset3=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]) & ((Yalt[:,0]==1) | (Yalt[:,4]==1)) )
#testset4=( ((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4])) & ((Yalt[:,2]==1) | (Yalt[:,4]==1)) )
#testset3=( ((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4])) & ((Yalt[:,0]==1) | (Yalt[:,4]==1)) )
#testset2=( ((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4])) & ((Yalt[:,1]==1) | (Yalt[:,4]==1)) )
#testset1=( ((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4])) & ((Yalt[:,3]==1) | (Yalt[:,4]==1)) )




#testset3=( ((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4])) & ((Yalt[:,0]==1) | (Yalt[:,1]==1)) )
#Pre = model.predict([Xchr[testset3,:,:],Xexp[testset3,:,:]])

#Yalt4=copy.deepcopy(Yalt)
#Yalt4=Yalt4[testset3,:]

#x1=Pre[:,[True,True,False,False,False]]
#x2=Yalt4[:,[True,True,False,False,False]]
#acc = sum([np.argmax(x1[i,:])==np.argmax(x2[i,:]) for i in range(np.shape(Yalt4)[0])])/np.shape(Yalt4)[0]

#Normalise the expression data
Xexp = copy.deepcopy(Xexpa)
Xexp[:,0,0] = ( Xexp[:,0,0] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0])
Xexp[:,0,1] = ( Xexp[:,0,1] - np.nanmean(Xexp[trainset,0,1]) ) / np.nanstd(Xexp[trainset,0,1])
#Xexp[:,0,3] = ( Xexp[:,0,3] - np.nanmean(Xexp[trainset,0,3]) ) / np.nanstd(Xexp[trainset,0,3])

Xchra = copy.deepcopy(Xchr)

#Normalise histone modification data
for k in range(0, np.shape(fil)[0]):
      Xchr[:,0:100,k] = ( Xchr[:,0:100,k] - np.nanmean(Xchr[:,0:100,k]) ) / np.nanstd(Xchr[:,0:100,k])


#Now load in the model
from keras.models import load_model
model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model1 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_neg.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model2 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_pos.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)

#classpred = model.predict([Xchr,Xexp], batch_size=1000)

Scores = np.zeros((3,1))
Scores[0,0] = model.evaluate([Xchr[trainset,:,:],Xexp[trainset,:,:]], Yalt[trainset,:], batch_size=32)[1]
Scores[1,0] = model.evaluate([Xchr[testset,:,:],Xexp[testset,:,:]], Yalt[testset,:], batch_size=32)[1]
Scores[2,0] = model.evaluate([Xchr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]

predictions = model.predict([Xchr,Xexp], batch_size=1000)
classpred = model.predict([Xchr,Xexp], batch_size=1000)


#NOW GENERATE SOME EXAMPLE PLOTS
#import matplotlib.pyplot as plt 
#from matplotlib import colors as mcolors
#plt.switch_backend('agg')
#from deepexplain.tensorflow import DeepExplain
#from keras import backend as K
#from sklearn.metrics import roc_curve, precision_recall_curve,auc

#del(Output2)
#del(Output)
#del(Xchra)
#Now do some heatmaps ...

#from deepexplain.tensorflow import DeepExplain
#from keras import backend as K




#Now do some others ...
#classpred = model1.predict([Xchr], batch_size=1000)

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
plt.switch_backend('agg')

#NoW PLOT THE SCATTER ploTs
uk, uinds = np.unique(Output2['f3'], return_index = True)
classpred1 = copy.deepcopy(classpred)
Yalt1 = copy.deepcopy(Yalt)
classpred1 = classpred1[uinds,:]
Yalt1 = Yalt1[uinds,:]
Xchr1 = copy.deepcopy(Xchr)
Xchr1 = Xchr1[uinds,:,:]
Xexp1 = copy.deepcopy(Xexp)
Xexp1 = Xexp1[uinds,:,:]
Xrawexp1 = copy.deepcopy(Xrawexp)
Xrawexp1 = Xrawexp1[uinds,:]
OnP = ((classpred1[:,0]>np.max(classpred1[:,1:5],axis=1)))
OnT = (Yalt1[:,0]==1)
OnTP = ((classpred1[:,0]>np.max(classpred1[:,1:5],axis=1)) & (Yalt1[:,0]==1))
OnFP = ((classpred1[:,0]>np.max(classpred1[:,1:5],axis=1)) & (Yalt1[:,0]==0))
OnTN = ((classpred1[:,0]<np.max(classpred1[:,1:5],axis=1)) & (Yalt1[:,0]==0))
OnFN = ((classpred1[:,0]<np.max(classpred1[:,1:5],axis=1)) & (Yalt1[:,0]==1))
OffT = (Yalt1[:,2]==1)
OffP = ((classpred1[:,2]>np.max(classpred1[:,[True,True,False,True,True]],axis=1) ))
OffTP = ((classpred1[:,2]>np.max(classpred1[:,[True,True,False,True,True]],axis=1)) & (Yalt1[:,2]==1))
OffFP = ((classpred1[:,2]>np.max(classpred1[:,[True,True,False,True,True]],axis=1)) & (Yalt1[:,2]==0))
OffTN = ((classpred1[:,2]<np.max(classpred1[:,[True,True,False,True,True]],axis=1)) & (Yalt1[:,2]==0))
OffFN = ((classpred1[:,2]<np.max(classpred1[:,[True,True,False,True,True]],axis=1)) & (Yalt1[:,2]==1))
RUT = (Yalt1[:,3]==1)
RUP = ((classpred1[:,3]>np.max(classpred1[:,[True,True,True,False,True]],axis=1) ))
RUTP = ((classpred1[:,3]>np.max(classpred1[:,[True,True,True,False,True]],axis=1)) & (Yalt1[:,3]==1))
RUFP = ((classpred1[:,3]>np.max(classpred1[:,[True,True,True,False,True]],axis=1)) & (Yalt1[:,3]==0))
RUTN = ((classpred1[:,3]<np.max(classpred1[:,[True,True,True,False,True]],axis=1)) & (Yalt1[:,3]==0))
RUFN = ((classpred1[:,3]<np.max(classpred1[:,[True,True,True,False,True]],axis=1)) & (Yalt1[:,3]==1))
RDT = (Yalt1[:,1]==1)
RDP = ((classpred1[:,1]>np.max(classpred1[:,[True,False,True,True,True]],axis=1) ))
RDTP = ((classpred1[:,1]>np.max(classpred1[:,[True,False,True,True,True]],axis=1)) & (Yalt1[:,1]==1))
RDFP = ((classpred1[:,1]>np.max(classpred1[:,[True,False,True,True,True]],axis=1)) & (Yalt1[:,1]==0))
RDTN = ((classpred1[:,1]<np.max(classpred1[:,[True,False,True,True,True]],axis=1)) & (Yalt1[:,1]==0))
RDFN = ((classpred1[:,1]<np.max(classpred1[:,[True,False,True,True,True]],axis=1)) & (Yalt1[:,1]==1))
OtherP = ((classpred1[:,4]>np.max(classpred1[:,[True,True,True,True,False]],axis=1) ))
OtherT = (Yalt1[:,4]==1)


cOtherT = np.zeros((np.shape(Xexp1[OtherT,:,:])[0],4))
cOtherT[:,0] = colors.to_rgba('lightgray')[0] #0.8274509803921568
cOtherT[:,1] = colors.to_rgba('lightgray')[1] #0.8274509803921568
cOtherT[:,2] = colors.to_rgba('lightgray')[2] #0.8274509803921568
cOtherT[:,3] = classpred1[OtherT,4]
cOffTP = np.zeros((np.shape(Xexp1[OffTP,:,:])[0],4))
cOffTP[:,0] = colors.to_rgba('royalblue')[0] #
cOffTP[:,1] = colors.to_rgba('royalblue')[1] #
cOffTP[:,2] = colors.to_rgba('royalblue')[2] #
cOffTP[:,3] = classpred1[OffTP,2]
cOnTP = np.zeros((np.shape(Xexp1[OnTP,:,:])[0],4))
cOnTP[:,0] = colors.to_rgba('orange')[0] #1.0
cOnTP[:,1] = colors.to_rgba('orange')[1] #0.5490196078431373
cOnTP[:,2] = colors.to_rgba('orange')[2] #0.5490196078431373
cOnTP[:,3] = classpred1[OnTP,0]
cOffFP = np.zeros((np.shape(Xexp1[OffFP,:,:])[0],4))
cOffFP[:,0] = colors.to_rgba('lightsteelblue')[0] #
cOffFP[:,1] = colors.to_rgba('lightsteelblue')[1] #
cOffFP[:,2] = colors.to_rgba('lightsteelblue')[2] #
cOffFP[:,3] = classpred1[OffFP,2]
cOnFP = np.zeros((np.shape(Xexp1[OnFP,:,:])[0],4))
cOnFP[:,0] = colors.to_rgba('bisque')[0] #1.0
cOnFP[:,1] = colors.to_rgba('bisque')[1] #0.5490196078431373
cOnFP[:,2] = colors.to_rgba('bisque')[2] #0.5490196078431373
cOnFP[:,3] = classpred1[OnFP,0]
cRDTP = np.zeros((np.shape(Xexp1[RDTP,:,:])[0],4))
cRDTP[:,0] = colors.to_rgba('bisque')[0]
cRDTP[:,1] = colors.to_rgba('bisque')[1]
cRDTP[:,2] = colors.to_rgba('bisque')[2]
cRDTP[:,3] = classpred1[RDTP,1]
cRUTP = np.zeros((np.shape(Xexp1[RUTP,:,:])[0],4))
cRUTP[:,0] = colors.to_rgba('lightsteelblue')[0]
cRUTP[:,1] = colors.to_rgba('lightsteelblue')[1]
cRUTP[:,2] = colors.to_rgba('lightsteelblue')[2]
cRUTP[:,3] = classpred1[RUTP,3]
cOffT = np.zeros((np.shape(Xexp1[OffT,:,:])[0],4))
cOffT[:,0] = colors.to_rgba('blue')[0] #
cOffT[:,1] = colors.to_rgba('blue')[1] #
cOffT[:,2] = colors.to_rgba('blue')[2] #
cOffT[:,3] = classpred1[OffT,2]
cOnT = np.zeros((np.shape(Xexp1[OnT,:,:])[0],4))
cOnT[:,0] = colors.to_rgba('red')[0] #1.0
cOnT[:,1] = colors.to_rgba('red')[1] #0.5490196078431373
cOnT[:,2] = colors.to_rgba('red')[2] #0.5490196078431373
cOnT[:,3] = classpred1[OnT,0]

#First plot everython coloured by p Xrawexp1
plt.figure()
plt.scatter((Xrawexp1[OtherT,5]+1),-np.log2(Xrawexp1[OtherT,6]), marker='o', color=cOtherT, linestyle='None')
plt.scatter((Xrawexp1[OnFP,5]+1),-np.log2(Xrawexp1[OnFP,6]), marker='o', color=cOnFP, linestyle='None')
plt.scatter((Xrawexp1[OnTP,5]+1),-np.log2(Xrawexp1[OnTP,6]), marker='o', color=cOnTP, linestyle='None')
plt.plot(np.linspace(-2,7,100),-np.log2(0.05)*np.ones((100)), 'k--')
plt.xlim((-8,8))
plt.ylim((-1,140))
plt.savefig('EndoFull/Endo_OnTPFP.pdf')
plt.draw()
plt.close()

#First plot everython coloured by p Xrawexp1
plt.figure()
plt.scatter((Xrawexp1[OtherT,5]+1),-np.log2(Xrawexp1[OtherT,6]), marker='o', color=cOtherT, linestyle='None')
plt.scatter((Xrawexp1[OffFP,5]+1),-np.log2(Xrawexp1[OffFP,6]), marker='o', color=cOffFP, linestyle='None')
plt.scatter((Xrawexp1[OffTP,5]+1),-np.log2(Xrawexp1[OffTP,6]), marker='o', color=cOffTP, linestyle='None')
plt.plot(np.linspace(-2,7,100),-np.log2(0.05)*np.ones((100)), 'k--')
plt.xlim((-8,8))
plt.ylim((-1,140))
plt.savefig('EndoFull/Endo_OffTPFP_chronly.pdf')
plt.draw()
plt.close()

#First plot everython coloured by p Xrawexp1
plt.figure()
plt.scatter((Xrawexp1[OtherT,5]+1),-np.log2(Xrawexp1[OtherT,6]), marker='o', color=cOtherT, linestyle='None')
plt.scatter((Xrawexp1[OffT,5]+1),-np.log2(Xrawexp1[OffT,6]), marker='o', color=cOffT, linestyle='None')
plt.scatter((Xrawexp1[OnT,5]+1),-np.log2(Xrawexp1[OnT,6]), marker='o', color=cOnT, linestyle='None')
plt.plot(np.linspace(-2,7,100),-np.log2(0.05)*np.ones((100)), 'k--')
plt.xlim((-8,8))
plt.ylim((-1,140))
plt.savefig('EndoFull/Endo_OnOffT.pdf')
plt.draw()
plt.close()

#Now do half models


#classpred = model.predict([Xchr,Xexp], batch_size=1000)
classpred = model1.predict([Xchr], batch_size=1000)

#from keras.models import load_model
#model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper.h5') #predictions = model.predict([Xchr,Xexp], batch_size=10$
#model1 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_neg.h5') #predictions = model.predict([Xchr,Xexp], batch_si$
#model2 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_pos.h5') #predictions = model.predict([Xchr,Xexp], batch_si$
#NOw SCATTER
uk, uinds = np.unique(Output2['f3'], return_index = True)
classpred1 = copy.deepcopy(classpred)
Yalt1 = copy.deepcopy(Yalt)
classpred1 = classpred1[uinds,:]
Yalt1 = Yalt1[uinds,:]
Xchr1 = copy.deepcopy(Xchr)
Xchr1 = Xchr1[uinds,:,:]
Xexp1 = copy.deepcopy(Xexp)
Xexp1 = Xexp1[uinds,:,:]
Xrawexp1 = copy.deepcopy(Xrawexp)
Xrawexp1 = Xrawexp1[uinds,:]
OnP = ((classpred1[:,0]>np.max(classpred1[:,1:5],axis=1)))
OnT = (Yalt1[:,0]==1)
OnTP = ((classpred1[:,0]>np.max(classpred1[:,1:5],axis=1)) & (Yalt1[:,0]==1))
OnFP = ((classpred1[:,0]>np.max(classpred1[:,1:5],axis=1)) & (Yalt1[:,0]==0))
OnTN = ((classpred1[:,0]<np.max(classpred1[:,1:5],axis=1)) & (Yalt1[:,0]==0))
OnFN = ((classpred1[:,0]<np.max(classpred1[:,1:5],axis=1)) & (Yalt1[:,0]==1))
OffT = (Yalt1[:,2]==1)
OffP = ((classpred1[:,2]>np.max(classpred1[:,[True,True,False,True,True]],axis=1) ))
OffTP = ((classpred1[:,2]>np.max(classpred1[:,[True,True,False,True,True]],axis=1)) & (Yalt1[:,2]==1))
OffFP = ((classpred1[:,2]>np.max(classpred1[:,[True,True,False,True,True]],axis=1)) & (Yalt1[:,2]==0))
OffTN = ((classpred1[:,2]<np.max(classpred1[:,[True,True,False,True,True]],axis=1)) & (Yalt1[:,2]==0))
OffFN = ((classpred1[:,2]<np.max(classpred1[:,[True,True,False,True,True]],axis=1)) & (Yalt1[:,2]==1))
RUT = (Yalt1[:,3]==1)
RUP = ((classpred1[:,3]>np.max(classpred1[:,[True,True,True,False,True]],axis=1) ))
RUTP = ((classpred1[:,3]>np.max(classpred1[:,[True,True,True,False,True]],axis=1)) & (Yalt1[:,3]==1))
RUFP = ((classpred1[:,3]>np.max(classpred1[:,[True,True,True,False,True]],axis=1)) & (Yalt1[:,3]==0))
RUTN = ((classpred1[:,3]<np.max(classpred1[:,[True,True,True,False,True]],axis=1)) & (Yalt1[:,3]==0))
RUFN = ((classpred1[:,3]<np.max(classpred1[:,[True,True,True,False,True]],axis=1)) & (Yalt1[:,3]==1))
RDT = (Yalt1[:,1]==1)
RDP = ((classpred1[:,1]>np.max(classpred1[:,[True,False,True,True,True]],axis=1) ))
RDTP = ((classpred1[:,1]>np.max(classpred1[:,[True,False,True,True,True]],axis=1)) & (Yalt1[:,1]==1))
RDFP = ((classpred1[:,1]>np.max(classpred1[:,[True,False,True,True,True]],axis=1)) & (Yalt1[:,1]==0))
RDTN = ((classpred1[:,1]<np.max(classpred1[:,[True,False,True,True,True]],axis=1)) & (Yalt1[:,1]==0))
RDFN = ((classpred1[:,1]<np.max(classpred1[:,[True,False,True,True,True]],axis=1)) & (Yalt1[:,1]==1))
OtherP = ((classpred1[:,4]>np.max(classpred1[:,[True,True,True,True,False]],axis=1) ))
OtherT = (Yalt1[:,4]==1)

cOtherT = np.zeros((np.shape(Xexp1[OtherT,:,:])[0],4))
cOtherT[:,0] = colors.to_rgba('lightgray')[0] #0.8274509803921568
cOtherT[:,1] = colors.to_rgba('lightgray')[1] #0.8274509803921568
cOtherT[:,2] = colors.to_rgba('lightgray')[2] #0.8274509803921568
cOtherT[:,3] = classpred1[OtherT,4]
cOffTP = np.zeros((np.shape(Xexp1[OffTP,:,:])[0],4))
cOffTP[:,0] = colors.to_rgba('royalblue')[0] #
cOffTP[:,1] = colors.to_rgba('royalblue')[1] #
cOffTP[:,2] = colors.to_rgba('royalblue')[2] #
cOffTP[:,3] = classpred1[OffTP,2]
cOnTP = np.zeros((np.shape(Xexp1[OnTP,:,:])[0],4))
cOnTP[:,0] = colors.to_rgba('orange')[0] #1.0
cOnTP[:,1] = colors.to_rgba('orange')[1] #0.5490196078431373
cOnTP[:,2] = colors.to_rgba('orange')[2] #0.5490196078431373
cOnTP[:,3] = classpred1[OnTP,0]
cOffFP = np.zeros((np.shape(Xexp1[OffFP,:,:])[0],4))
cOffFP[:,0] = colors.to_rgba('lightsteelblue')[0] #
cOffFP[:,1] = colors.to_rgba('lightsteelblue')[1] #
cOffFP[:,2] = colors.to_rgba('lightsteelblue')[2] #
cOffFP[:,3] = classpred1[OffFP,2]
cOnFP = np.zeros((np.shape(Xexp1[OnFP,:,:])[0],4))
cOnFP[:,0] = colors.to_rgba('bisque')[0] #1.0
cOnFP[:,1] = colors.to_rgba('bisque')[1] #0.5490196078431373
cOnFP[:,2] = colors.to_rgba('bisque')[2] #0.5490196078431373
cOnFP[:,3] = classpred1[OnFP,0]
cRDTP = np.zeros((np.shape(Xexp1[RDTP,:,:])[0],4))
cRDTP[:,0] = colors.to_rgba('bisque')[0]
cRDTP[:,1] = colors.to_rgba('bisque')[1]
cRDTP[:,2] = colors.to_rgba('bisque')[2]
cRDTP[:,3] = classpred1[RDTP,1]
cRUTP = np.zeros((np.shape(Xexp1[RUTP,:,:])[0],4))
cRUTP[:,0] = colors.to_rgba('lightsteelblue')[0]
cRUTP[:,1] = colors.to_rgba('lightsteelblue')[1]
cRUTP[:,2] = colors.to_rgba('lightsteelblue')[2]
cRUTP[:,3] = classpred1[RUTP,3]
cOffT = np.zeros((np.shape(Xexp1[OffT,:,:])[0],4))
cOffT[:,0] = colors.to_rgba('blue')[0] #
cOffT[:,1] = colors.to_rgba('blue')[1] #
cOffT[:,2] = colors.to_rgba('blue')[2] #
cOffT[:,3] = classpred1[OffT,2]
cOnT = np.zeros((np.shape(Xexp1[OnT,:,:])[0],4))
cOnT[:,0] = colors.to_rgba('red')[0] #1.0
cOnT[:,1] = colors.to_rgba('red')[1] #0.5490196078431373
cOnT[:,2] = colors.to_rgba('red')[2] #0.5490196078431373
cOnT[:,3] = classpred1[OnT,0]

#First plot everython coloured by p Xrawexp1
plt.figure()
plt.scatter((Xrawexp1[OtherT,5]+1),-np.log2(Xrawexp1[OtherT,6]), marker='o', color=cOtherT, linestyle='None')
plt.scatter((Xrawexp1[OnFP,5]+1),-np.log2(Xrawexp1[OnFP,6]), marker='o', color=cOnFP, linestyle='None')
plt.scatter((Xrawexp1[OnTP,5]+1),-np.log2(Xrawexp1[OnTP,6]), marker='o', color=cOnTP, linestyle='None')
plt.plot(np.linspace(-2,7,100),-np.log2(0.05)*np.ones((100)), 'k--')
plt.xlim((-8,8))
plt.ylim((-1,140))
plt.savefig('EndoChr/Endo_OnTPFP_chronly.pdf')
plt.draw()
plt.close()

#First plot everython coloured by p Xrawexp1
plt.figure()
plt.scatter((Xrawexp1[OtherT,5]+1),-np.log2(Xrawexp1[OtherT,6]), marker='o', color=cOtherT, linestyle='None')
plt.scatter((Xrawexp1[OffFP,5]+1),-np.log2(Xrawexp1[OffFP,6]), marker='o', color=cOffFP, linestyle='None')
plt.scatter((Xrawexp1[OffTP,5]+1),-np.log2(Xrawexp1[OffTP,6]), marker='o', color=cOffTP, linestyle='None')
plt.plot(np.linspace(-2,7,100),-np.log2(0.05)*np.ones((100)), 'k--')
plt.xlim((-8,8))
plt.ylim((-1,140))
plt.savefig('EndoChr/Endo_OffTPFP_chronly.pdf')
plt.draw()
plt.close()

#First plot everython coloured by p Xrawexp1
plt.figure()
plt.scatter((Xrawexp1[OtherT,5]+1),-np.log2(Xrawexp1[OtherT,6]), marker='o', color=cOtherT, linestyle='None')
plt.scatter((Xrawexp1[OffT,5]+1),-np.log2(Xrawexp1[OffT,6]), marker='o', color=cOffT, linestyle='None')
plt.scatter((Xrawexp1[OnT,5]+1),-np.log2(Xrawexp1[OnT,6]), marker='o', color=cOnT, linestyle='None')
plt.plot(np.linspace(-2,7,100),-np.log2(0.05)*np.ones((100)), 'k--')
plt.xlim((-8,8))
plt.ylim((-1,140))
plt.savefig('EndoChr/Endo_OnOffT.pdf')
plt.draw()
plt.close()





model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model2 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_neg.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model3 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_pos.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model4 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model5 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_neg_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model6 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_pos_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)

#Endoderm/Model_deeper.h5		Endoderm/Model_deeper_minimalset.h5	Endoderm/Model_deeper_neg.h5		Endoderm/Model_deeper_pos.h5
Scores = np.zeros((3,1))
Scores[0,0] = model.evaluate([Xchr[trainset,:,:],Xexp[trainset,:,:]], Yalt[trainset,:], batch_size=32)[1]
Scores[1,0] = model.evaluate([Xchr[testset,:,:],Xexp[testset,:,:]], Yalt[testset,:], batch_size=32)[1]
Scores[2,0] = model.evaluate([Xchr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
classpred = model.predict([Xchr,Xexp], batch_size=1000)
classpred2 = model2.predict([Xchr], batch_size=1000)
classpred3 = model3.predict([Xexp], batch_size=1000)
classpred4 = model.predict([Xchr,Xexp], batch_size=1000)
classpred5 = model2.predict([Xchr], batch_size=1000)
classpred6 = model3.predict([Xexp], batch_size=1000)

#NOW GENERATE SOME EXAMPLE PLOTS
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
plt.switch_backend('agg')
#classpred1 = model.predict([Xchr[valset,:,:],,Xexp[valset,:,:],], batch_size=1000)

from sklearn.metrics import roc_curve, precision_recall_curve,auc
fpr_keras, tpr_keras, thresholds_keras = roc_curve(Yalt.ravel(), model.predict([Xchr,Xexp], batch_size=1000).ravel())
pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt.ravel(), model.predict([Xchr,Xexp], batch_size=1000).ravel())
fpr_kerast, tpr_kerast, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_kerast, re_kerast, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
fpr_keras1, tpr_keras1, thresholds_keras = roc_curve(Yalt[valset,:].ravel(), model2.predict(Xchr[valset,:,:], batch_size=1000).ravel())
pr_keras1, re_keras1, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model2.predict(Xchr[valset,:,:], batch_size=1000).ravel())
fpr_keras2, tpr_keras2, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model3.predict(Xexp[valset,:,:], batch_size=1000).ravel())
pr_keras2, re_keras2, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model3.predict(Xexp[valset,:,:], batch_size=1000).ravel())
fpr_kerast4, tpr_kerast4, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model4.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_kerast4, re_kerast4, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model4.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
fpr_keras5, tpr_keras5, thresholds_keras = roc_curve(Yalt[valset,:].ravel(), model5.predict(Xchr[valset,:,:], batch_size=1000).ravel())
pr_keras5, re_keras5, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model5.predict(Xchr[valset,:,:], batch_size=1000).ravel())
fpr_keras6, tpr_keras6, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model6.predict(Xexp[valset,:,:], batch_size=1000).ravel())
pr_keras6, re_keras6, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model6.predict(Xexp[valset,:,:], batch_size=1000).ravel())

plt.figure()
plt.plot(re_kerast, pr_kerast, label='Combined') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(re_keras1, pr_keras1, label='Histone') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(re_keras2, pr_keras2, label='Expression') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(re_kerast4, pr_kerast4, label='Combined (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(re_keras5, pr_keras5, label='Histone (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(re_keras6, pr_keras6, label='Expression (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
#plt.xlim([0.0, 1.0])
#plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig('EndoFull/PR.pdf')
plt.draw()

plt.figure()
plt.plot(fpr_kerast, tpr_kerast, label='Combined') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(fpr_keras1, tpr_keras1, label='Histone') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(fpr_keras2, tpr_keras2, label='Expression') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(fpr_kerast4, tpr_kerast4, label='Combined (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(fpr_keras5, tpr_keras5, label='Histone (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(fpr_keras6, tpr_keras6, label='Expression (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot([0, 1], [0, 1], 'k--', label='Random')
plt.plot([0,0, 1], [0, 1, 1], 'k-')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig('EndoFull/AUC.pdf')
plt.draw()









