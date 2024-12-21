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

np.save('EndoFull/Output2.npy', Output2)


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

np.save('EndoFull/Yalt.npy', Yalt)


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


#np.save('EndoFull/Yalt.npy', Xex)

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

#Normalise the expression data
Xexp = copy.deepcopy(Xexpa)
Xexp[:,0,0] = ( Xexp[:,0,0] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0])
Xexp[:,0,1] = ( Xexp[:,0,1] - np.nanmean(Xexp[trainset,0,1]) ) / np.nanstd(Xexp[trainset,0,1])
#Xexp[:,0,3] = ( Xexp[:,0,3] - np.nanmean(Xexp[trainset,0,3]) ) / np.nanstd(Xexp[trainset,0,3])

Xchra = copy.deepcopy(Xchr)

#Normalise histone modification data
for k in range(0, np.shape(fil)[0]):
      Xchr[:,0:600,k] = ( Xchr[:,0:600,k] - np.nanmean(Xchr[trainset,0:600,k]) ) / np.nanstd(Xchr[trainset,0:600,k])


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

np.save('EndoFull/OnTP.npy', OnTP)
np.save('EndoFull/OnFP.npy', OnFP)
np.save('EndoFull/OnTN.npy', OnTN)
np.save('EndoFull/OnFN.npy', OnFN)
np.save('EndoFull/OffTP.npy', OffTP)
np.save('EndoFull/OffFP.npy', OffFP)
np.save('EndoFull/OffTN.npy', OffTN)
np.save('EndoFull/OffFN.npy', OffFN)
np.save('EndoFull/RUTP.npy', RUTP)
np.save('EndoFull/RUFP.npy', RUFP)
np.save('EndoFull/RUTN.npy', RUTN)
np.save('EndoFull/RUFN.npy', RUFN)
np.save('EndoFull/RDTP.npy', RDTP)
np.save('EndoFull/RDFP.npy', RDFP)
np.save('EndoFull/RDTN.npy', RDTN)
np.save('EndoFull/RDFN.npy', RDFN)

xs1_1 = Xchr[OnTP,:,:]
xs2_1 = Xexp[OnTP,:,:]
ys_1 = Yalt[OnTP,:]

xs1_2 = Xchr[OffTP,:,:]
xs2_2 = Xexp[OffTP,:,:]
ys_2 = Yalt[OffTP,:]

xs1_3 = Xchr[RUTP,:,:]
xs2_3 = Xexp[RUTP,:,:]
ys_3 = Yalt[RUTP,:]

xs1_4 = Xchr[RDTP,:,:]
xs2_4 = Xexp[RDTP,:,:]
ys_4 = Yalt[RDTP,:]

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions4 = de.explain('grad*input', target_tensor*np.array([1., 0., 0., 0., 0.]), input_tensors, [xs1_1,xs2_1])
    attributions4B = de.explain('grad*input', target_tensor*np.array([0., 1., 0., 0., 0.]), input_tensors, [xs1_1,xs2_1])
    attributions4C = de.explain('grad*input', target_tensor*np.array([0., 0., 1., 0., 0.]), input_tensors, [xs1_1,xs2_1])
    attributions4D = de.explain('grad*input', target_tensor*np.array([0., 0., 0., 1., 0.]), input_tensors, [xs1_1,xs2_1])
    attributions4E = de.explain('grad*input', target_tensor*np.array([0., 0., 0., 0., 1.]), input_tensors, [xs1_1,xs2_1])

np.save('EndoFull/On_attributions_OnTP.npy', attributions4[0])
np.save('EndoFull/Down_attributions_OnTP.npy', attributions4B[0])

np.save('EndoFull/Off_attributions_OnTP.npy', attributions4C[0])
np.save('EndoFull/Up_attributions_OnTP.npy', attributions4D[0])
np.save('EndoFull/Other_attributions_OnTP.npy', attributions4E[0])


np.save('EndoFull/On_attributions_OnTP_1.npy', attributions4[1])
np.save('EndoFull/Down_attributions_OnTP_1.npy', attributions4B[1])
np.save('EndoFull/Off_attributions_OnTP_1.npy', attributions4C[1])
np.save('EndoFull/Up_attributions_OnTP_1.npy', attributions4D[1])
np.save('EndoFull/Other_attributions_OnTP_1.npy', attributions4E[1])

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions5 = de.explain('grad*input', target_tensor*np.array([0, 0., 1., 0., 0.]), input_tensors, [xs1_2,xs2_2])
    attributions5B = de.explain('grad*input', target_tensor*np.array([0, 0., 0., 1., 0.]), input_tensors, [xs1_2,xs2_2])
    attributions5C = de.explain('grad*input', target_tensor*np.array([1, 0., 0., 0., 0.]), input_tensors, [xs1_2,xs2_2])
    attributions5D = de.explain('grad*input', target_tensor*np.array([0, 1., 0., 0., 0.]), input_tensors, [xs1_2,xs2_2])
    attributions5E = de.explain('grad*input', target_tensor*np.array([0, 0., 0., 0., 1.]), input_tensors, [xs1_2,xs2_2])

np.save('EndoFull/Off_attributions_OffTP.npy', attributions5[0])
np.save('EndoFull/Up_attributions_OffTP.npy', attributions5B[0])
np.save('EndoFull/On_attributions_OffTP.npy', attributions5C[0])
np.save('EndoFull/Down_attributions_OffTP.npy', attributions5D[0])
np.save('EndoFull/Other_attributions_OffTP.npy', attributions5E[0])

np.save('EndoFull/Off_attributions_OffTP_1.npy', attributions5[1])
np.save('EndoFull/Up_attributions_OffTP_1.npy', attributions5B[1])
np.save('EndoFull/On_attributions_OffTP_1.npy', attributions5C[1])
np.save('EndoFull/Down_attributions_OffTP_1.npy', attributions5D[1])
np.save('EndoFull/Other_attributions_OffTP_1.npy', attributions5E[1])


with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions4 = de.explain('grad*input', target_tensor*np.array([0, 0., 0., 1., 0.]), input_tensors, [xs1_3,xs2_3])
    attributions4B = de.explain('grad*input', target_tensor*np.array([0, 0., 1., 0., 0.]), input_tensors, [xs1_3,xs2_3])
    attributions4C = de.explain('grad*input', target_tensor*np.array([0., 0., 1., 0., 0.]), input_tensors, [xs1_1,xs2_1])
    attributions4D = de.explain('grad*input', target_tensor*np.array([0., 0., 0., 1., 0.]), input_tensors, [xs1_1,xs2_1])
    attributions4E = de.explain('grad*input', target_tensor*np.array([0., 0., 0., 0., 1.]), input_tensors, [xs1_1,xs2_1])


np.save('EndoFull/Up_attributions_RUTP.npy', attributions4[0])
np.save('EndoFull/Off_attributions_RUTP.npy', attributions4B[0])
np.save('EndoFull/Up_attributions_RUTP_1.npy', attributions4[1])
np.save('EndoFull/Off_attributions_RUTP_1.npy', attributions4B[1])

np.save('EndoFull/Off_attributions_OnTP.npy', attributions4C[0])
np.save('EndoFull/Up_attributions_OnTP.npy', attributions4D[0])
np.save('EndoFull/Other_attributions_OnTP.npy', attributions4E[0])
np.save('EndoFull/Off_attributions_OnTP_1.npy', attributions4C[1])
np.save('EndoFull/Up_attributions_OnTP_1.npy', attributions4D[1])
np.save('EndoFull/Other_attributions_OnTP_1.npy', attributions4E[1])



with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions5 = de.explain('grad*input', target_tensor*np.array([0, 1., 0., 0., 0.]), input_tensors, [xs1_4,xs2_4])
    attributions5B = de.explain('grad*input', target_tensor*np.array([1, 0., 0., 0., 0.]), input_tensors, [xs1_4,xs2_4])
    attributions5C = de.explain('grad*input', target_tensor*np.array([1, 0., 0., 0., 0.]), input_tensors, [xs1_2,xs2_2])
    attributions5D = de.explain('grad*input', target_tensor*np.array([0, 1., 0., 0., 0.]), input_tensors, [xs1_2,xs2_2])
    attributions5E = de.explain('grad*input', target_tensor*np.array([0, 0., 0., 0., 1.]), input_tensors, [xs1_2,xs2_2])


np.save('EndoFull/Down_attributions_RDTP.npy', attributions5[0])
np.save('EndoFull/On_attributions_RDTP.npy', attributions5B[0])
np.save('EndoFull/Down_attributions_RDTP_1.npy', attributions5[1])
np.save('EndoFull/On_attributions_RDTP_1.npy', attributions5B[1])
np.save('EndoFull/On_attributions_OffTP.npy', attributions5C[0])
np.save('EndoFull/Down_attributions_OffTP.npy', attributions5D[0])
np.save('EndoFull/Other_attributions_OffTP.npy', attributions5E[0])
np.save('EndoFull/On_attributions_OffTP_1.npy', attributions5C[1])
np.save('EndoFull/Down_attributions_OffTP_1.npy', attributions5D[1])
np.save('EndoFull/Other_attributions_OffTP_1.npy', attributions5E[1])


OnTP = np.random.permutation(OnTP)
OffTP = np.random.permutation(OffTP)
RUTP = np.random.permutation(RUTP)
RDTP = np.random.permutation(RDTP)

np.save('EndoFull/OnTP_shuffle.npy', OnTP)
np.save('EndoFull/OffTP_shuffle.npy', OffTP)
np.save('EndoFull/RUTP_shuffle.npy', RUTP)
np.save('EndoFull/RDTP_shuffle.npy', RDTP)


xs1_1 = Xchr[OnTP,:,:]
xs2_1 = Xexp[OnTP,:,:]
ys_1 = Yalt[OnTP,:]

xs1_2 = Xchr[OffTP,:,:]
xs2_2 = Xexp[OffTP,:,:]
ys_2 = Yalt[OffTP,:]

xs1_3 = Xchr[RUTP,:,:]
xs2_3 = Xexp[RUTP,:,:]
ys_3 = Yalt[RUTP,:]

xs1_4 = Xchr[RDTP,:,:]
xs2_4 = Xexp[RDTP,:,:]
ys_4 = Yalt[RDTP,:]

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions4 = de.explain('grad*input', target_tensor*np.array([1., 0., 0., 0., 0.]), input_tensors, [xs1_1,xs2_1])
    attributions4B = de.explain('grad*input', target_tensor*np.array([0., 1., 0., 0., 0.]), input_tensors, [xs1_1,xs2_1])

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions5 = de.explain('grad*input', target_tensor*np.array([0, 0., 1., 0., 0.]), input_tensors, [xs1_2,xs2_2])
    attributions5B = de.explain('grad*input', target_tensor*np.array([0, 0., 0., 1., 0.]), input_tensors, [xs1_2,xs2_2])

np.save('EndoFull/On_attributions_OnTP_perm.npy', attributions4[0])
np.save('EndoFull/Down_attributions_OnTP_perm.npy', attributions4B[0])
np.save('EndoFull/Off_attributions_OffTP_perm.npy', attributions5[0])
np.save('EndoFull/Up_attributions_OffTP_perm.npy', attributions5B[0])

np.save('EndoFull/On_attributions_OnTP_1_perm.npy', attributions4[1])
np.save('EndoFull/Down_attributions_OnTP_1_perm.npy', attributions4B[1])
np.save('EndoFull/Off_attributions_OffTP_1_perm.npy', attributions5[1])
np.save('EndoFull/Up_attributions_OffTP_1_perm.npy', attributions5B[1])


with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions4 = de.explain('grad*input', target_tensor*np.array([0, 0., 0., 1., 0.]), input_tensors, [xs1_3,xs2_3])
    attributions4B = de.explain('grad*input', target_tensor*np.array([0, 0., 1., 0., 0.]), input_tensors, [xs1_3,xs2_3])

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions5 = de.explain('grad*input', target_tensor*np.array([0, 1., 0., 0., 0.]), input_tensors, [xs1_4,xs2_4])
    attributions5B = de.explain('grad*input', target_tensor*np.array([1, 0., 0., 0., 0.]), input_tensors, [xs1_4,xs2_4])

np.save('EndoFull/Up_attributions_RUTP_perm.npy', attributions4[0])
np.save('EndoFull/Off_attributions_RUTP_perm.npy', attributions4B[0])
np.save('EndoFull/Down_attributions_RDTP_perm.npy', attributions5[0])
np.save('EndoFull/On_attributions_RDTP_perm.npy', attributions5B[0])


np.save('EndoFull/Up_attributions_RUTP_1_perm.npy', attributions4[1])
np.save('EndoFull/Off_attributions_RUTP_1_perm.npy', attributions4B[1])
np.save('EndoFull/Down_attributions_RDTP_1_perm.npy', attributions5[1])
np.save('EndoFull/On_attributions_RDTP_1_perm.npy', attributions5B[1])

from sklearn.ensemble import RandomForestClassifier


A =np.squeeze(np.reshape(Xchr[trainset,:,:],(np.shape(Xchr[trainset,:,:])[0],   np.shape(Xchr[trainset,:,:])[1] *np.shape(Xchr[trainset,:,:])[2])  ))
B =np.squeeze(Xexp[trainset,:,:] )
C = np.concatenate((A,B),axis=1)

clf = RandomForestClassifier(n_estimators=200, max_depth=2,random_state=0)
clf.fit(C , np.argmax(Yalt[trainset,:], axis=1) )
A1 =np.squeeze(np.reshape(Xchr[valset,:,:],(np.shape(Xchr[valset,:,:])[0],   np.shape(Xchr[valset,:,:])[1] *np.shape(Xchr[valset,:,:])[2])  ))
B1 =np.squeeze(Xexp[valset,:,:] )
C1 = np.concatenate((A1,B1),axis=1)
rfpred = clf.predict_proba(C1)
np.save('EndoFull/rfpreds1.npy',rfpred)
fpr_keras, tpr_keras, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), rfpred.ravel())
pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), rfpred.ravel())
np.save('EndoFull/rf_fpr_kearss1.npy', fpr_keras)
np.save('EndoFull/rf_tpr_kearss1.npy', tpr_keras)
np.save('EndoFull/rf_pr_kearss1.npy', pr_keras)
np.save('EndoFull/rf_re_kearss1.npy', re_keras)

P0 = rfpred
clf = RandomForestClassifier(n_estimators=200, max_depth=2,random_state=1)
clf.fit(C , np.argmax(Yalt[trainset,:], axis=1) )
A1 =np.squeeze(np.reshape(Xchr[valset,:,:],(np.shape(Xchr[valset,:,:])[0],   np.shape(Xchr[valset,:,:])[1] *np.shape(Xchr[valset,:,:])[2])  ))
B1 =np.squeeze(Xexp[valset,:,:] )
C1 = np.concatenate((A1,B1),axis=1)
rfpred = clf.predict_proba(C1)
np.save('EndoFull/rfpreds2.npy',rfpred)
fpr_keras, tpr_keras, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), rfpred.ravel())
pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), rfpred.ravel())
np.save('EndoFull/rf_fpr_kearss2.npy', fpr_keras)
np.save('EndoFull/rf_tpr_kearss2.npy', tpr_keras)
np.save('EndoFull/rf_pr_kearss2.npy', pr_keras)
np.save('EndoFull/rf_re_kearss2.npy', re_keras)
P1 = rfpred
Y0 = Yalt[valset,:]
Y1 = Yalt[valset,:]
fpr_keras0_1, tpr_keras0_1, thresholds_kerast = roc_curve(Y0[:,[True,False,False,False,False]], P0[:,[True,False,False,False,False]])
fpr_keras0_2, tpr_keras0_2, thresholds_kerast = roc_curve(Y0[:,[False,True,False,False,False]], P0[:,[False,True,False,False,False]])
fpr_keras0_3, tpr_keras0_3, thresholds_kerast = roc_curve(Y0[:,[False,False,True,False,False]], P0[:,[False,False,True,False,False]])
fpr_keras0_4, tpr_keras0_4, thresholds_kerast = roc_curve(Y0[:,[False,False,False,True,False]], P0[:,[False,False,False,True,False]])
fpr_keras0_5, tpr_keras0_5, thresholds_kerast = roc_curve(Y0[:,[False,False,False,False,True]], P0[:,[False,False,False,False,True]])
fpr_keras1_1, tpr_keras1_1, thresholds_kerast = roc_curve(Y1[:,[True,False,False,False,False]], P1[:,[True,False,False,False,False]])
fpr_keras1_2, tpr_keras1_2, thresholds_kerast = roc_curve(Y1[:,[False,True,False,False,False]], P1[:,[False,True,False,False,False]])
fpr_keras1_3, tpr_keras1_3, thresholds_kerast = roc_curve(Y1[:,[False,False,True,False,False]], P1[:,[False,False,True,False,False]])
fpr_keras1_4, tpr_keras1_4, thresholds_kerast = roc_curve(Y1[:,[False,False,False,True,False]], P1[:,[False,False,False,True,False]])
fpr_keras1_5, tpr_keras1_5, thresholds_kerast = roc_curve(Y1[:,[False,False,False,False,True]], P1[:,[False,False,False,False,True]])
np.save('EndoFull/rf_tpr_kears0_1.npy', tpr_keras0_1)
np.save('EndoFull/rf_tpr_kears0_2.npy', tpr_keras0_2)
np.save('EndoFull/rf_tpr_kears0_3.npy', tpr_keras0_3)
np.save('EndoFull/rf_tpr_kears0_4.npy', tpr_keras0_4)
np.save('EndoFull/rf_tpr_kears0_5.npy', tpr_keras0_5)
np.save('EndoFull/rf_tpr_kears1_1.npy', tpr_keras1_1)
np.save('EndoFull/rf_tpr_kears1_2.npy', tpr_keras1_2)
np.save('EndoFull/rf_tpr_kears1_3.npy', tpr_keras1_3)
np.save('EndoFull/rf_tpr_kears1_4.npy', tpr_keras1_4)
np.save('EndoFull/rf_tpr_kears1_5.npy', tpr_keras1_5)
np.save('EndoFull/rf_fpr_kears0_1.npy', fpr_keras0_1)
np.save('EndoFull/rf_fpr_kears0_2.npy', fpr_keras0_2)
np.save('EndoFull/rf_fpr_kears0_3.npy', fpr_keras0_3)
np.save('EndoFull/rf_fpr_kears0_4.npy', fpr_keras0_4)
np.save('EndoFull/rf_fpr_kears0_5.npy', fpr_keras0_5)
np.save('EndoFull/rf_fpr_kears1_1.npy', fpr_keras1_1)
np.save('EndoFull/rf_fpr_kears1_2.npy', fpr_keras1_2)
np.save('EndoFull/rf_fpr_kears1_3.npy', fpr_keras1_3)
np.save('EndoFull/rf_fpr_kears1_4.npy', fpr_keras1_4)
np.save('EndoFull/rf_fpr_kears1_5.npy', fpr_keras1_5)
pr_keras0_1, re_keras0_1, thresholds_kerast = precision_recall_curve(Y0[:,[True,False,False,False,False]], P0[:,[True,False,False,False,False]])
pr_keras0_2, re_keras0_2, thresholds_kerast = precision_recall_curve(Y0[:,[False,True,False,False,False]], P0[:,[False,True,False,False,False]])
pr_keras0_3, re_keras0_3, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,True,False,False]], P0[:,[False,False,True,False,False]])
pr_keras0_4, re_keras0_4, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,False,True,False]], P0[:,[False,False,False,True,False]])
pr_keras0_5, re_keras0_5, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,False,False,True]], P0[:,[False,False,False,False,True]])
pr_keras1_1, re_keras1_1, thresholds_kerast = precision_recall_curve(Y1[:,[True,False,False,False,False]], P1[:,[True,False,False,False,False]])
pr_keras1_2, re_keras1_2, thresholds_kerast = precision_recall_curve(Y1[:,[False,True,False,False,False]], P1[:,[False,True,False,False,False]])
pr_keras1_3, re_keras1_3, thresholds_kerast = precision_recall_curve(Y1[:,[False,False,True,False,False]], P1[:,[False,False,True,False,False]])
pr_keras1_4, re_keras1_4, thresholds_kerast = precision_recall_curve(Y1[:,[False,False,False,True,False]], P1[:,[False,False,False,True,False]])
pr_keras1_5, re_keras1_5, thresholds_kerast = precision_recall_curve(Y1[:,[False,False,False,False,True]], P1[:,[False,False,False,False,True]])
np.save('EndoFull/rf_pr_kears0_1.npy', pr_keras0_1)
np.save('EndoFull/rf_pr_kears0_2.npy', pr_keras0_2)
np.save('EndoFull/rf_pr_kears0_3.npy', pr_keras0_3)
np.save('EndoFull/rf_pr_kears0_4.npy', pr_keras0_4)
np.save('EndoFull/rf_pr_kears0_5.npy', pr_keras0_5)
np.save('EndoFull/rf_pr_kears1_1.npy', pr_keras1_1)
np.save('EndoFull/rf_pr_kears1_2.npy', pr_keras1_2)
np.save('EndoFull/rf_pr_kears1_3.npy', pr_keras1_3)
np.save('EndoFull/rf_pr_kears1_4.npy', pr_keras1_4)
np.save('EndoFull/rf_pr_kears1_5.npy', pr_keras1_5)
np.save('EndoFull/rf_re_kears0_1.npy', re_keras0_1)
np.save('EndoFull/rf_re_kears0_2.npy', re_keras0_2)
np.save('EndoFull/rf_re_kears0_3.npy', re_keras0_3)
np.save('EndoFull/rf_re_kears0_4.npy', re_keras0_4)
np.save('EndoFull/rf_re_kears0_5.npy', re_keras0_5)
np.save('EndoFull/rf_re_kears1_1.npy', re_keras1_1)
np.save('EndoFull/rf_re_kears1_2.npy', re_keras1_2)
np.save('EndoFull/rf_re_kears1_3.npy', re_keras1_3)
np.save('EndoFull/rf_re_kears1_4.npy', re_keras1_4)
np.save('EndoFull/rf_re_kears1_5.npy', re_keras1_5)


A =np.squeeze(np.reshape(Xchr[trainset,:,:],(np.shape(Xchr[trainset,:,:])[0],   np.shape(Xchr[trainset,:,:])[1] *np.shape(Xchr[trainset,:,:])[2])  ))
B =np.squeeze(Xexp[trainset,:,:] )
C = np.concatenate((A,B),axis=1)

clf = RandomForestClassifier(n_estimators=200, max_depth=2,random_state=0)
clf.fit(B , np.argmax(Yalt[trainset,:], axis=1) )
A1 =np.squeeze(np.reshape(Xchr[valset,:,:],(np.shape(Xchr[valset,:,:])[0],   np.shape(Xchr[valset,:,:])[1] *np.shape(Xchr[valset,:,:])[2])  ))
B1 =np.squeeze(Xexp[valset,:,:] )
C1 = np.concatenate((A1,B1),axis=1)
rfpred = clf.predict_proba(B1)
np.save('EndoFull/rfpreds1_neg.npy',rfpred)
fpr_keras, tpr_keras, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), rfpred.ravel())
pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), rfpred.ravel())
np.save('EndoFull/rf_fpr_kearss1_neg.npy', fpr_keras)
np.save('EndoFull/rf_tpr_kearss1_neg.npy', tpr_keras)
np.save('EndoFull/rf_pr_kearss1_neg.npy', pr_keras)
np.save('EndoFull/rf_re_kearss1_neg.npy', re_keras)

P0 = rfpred
clf = RandomForestClassifier(n_estimators=200, max_depth=2,random_state=1)
clf.fit(B , np.argmax(Yalt[trainset,:], axis=1) )
A1 =np.squeeze(np.reshape(Xchr[valset,:,:],(np.shape(Xchr[valset,:,:])[0],   np.shape(Xchr[valset,:,:])[1] *np.shape(Xchr[valset,:,:])[2])  ))
B1 =np.squeeze(Xexp[valset,:,:] )
C1 = np.concatenate((A1,B1),axis=1)
rfpred = clf.predict_proba(B1)
np.save('EndoFull/rfpreds2_neg.npy',rfpred)
fpr_keras, tpr_keras, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), rfpred.ravel())
pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), rfpred.ravel())
np.save('EndoFull/rf_fpr_kearss2_neg.npy', fpr_keras)
np.save('EndoFull/rf_tpr_kearss2_neg.npy', tpr_keras)
np.save('EndoFull/rf_pr_kearss2_neg.npy', pr_keras)
np.save('EndoFull/rf_re_kearss2_neg.npy', re_keras)
P1 = rfpred
Y0 = Yalt[valset,:]
Y1 = Yalt[valset,:]
fpr_keras0_1, tpr_keras0_1, thresholds_kerast = roc_curve(Y0[:,[True,False,False,False,False]], P0[:,[True,False,False,False,False]])
fpr_keras0_2, tpr_keras0_2, thresholds_kerast = roc_curve(Y0[:,[False,True,False,False,False]], P0[:,[False,True,False,False,False]])
fpr_keras0_3, tpr_keras0_3, thresholds_kerast = roc_curve(Y0[:,[False,False,True,False,False]], P0[:,[False,False,True,False,False]])
fpr_keras0_4, tpr_keras0_4, thresholds_kerast = roc_curve(Y0[:,[False,False,False,True,False]], P0[:,[False,False,False,True,False]])
fpr_keras0_5, tpr_keras0_5, thresholds_kerast = roc_curve(Y0[:,[False,False,False,False,True]], P0[:,[False,False,False,False,True]])
fpr_keras1_1, tpr_keras1_1, thresholds_kerast = roc_curve(Y1[:,[True,False,False,False,False]], P1[:,[True,False,False,False,False]])
fpr_keras1_2, tpr_keras1_2, thresholds_kerast = roc_curve(Y1[:,[False,True,False,False,False]], P1[:,[False,True,False,False,False]])
fpr_keras1_3, tpr_keras1_3, thresholds_kerast = roc_curve(Y1[:,[False,False,True,False,False]], P1[:,[False,False,True,False,False]])
fpr_keras1_4, tpr_keras1_4, thresholds_kerast = roc_curve(Y1[:,[False,False,False,True,False]], P1[:,[False,False,False,True,False]])
fpr_keras1_5, tpr_keras1_5, thresholds_kerast = roc_curve(Y1[:,[False,False,False,False,True]], P1[:,[False,False,False,False,True]])

fpr_keras0_1, tpr_keras0_1, thresholds_kerast = roc_curve(Y0[:,[True,False,False,False,False]], P0[:,[True,False,False,False,False]])
fpr_keras0_2, tpr_keras0_2, thresholds_kerast = roc_curve(Y0[:,[False,True,False,False,False]], P0[:,[False,True,False,False,False]])
fpr_keras0_3, tpr_keras0_3, thresholds_kerast = roc_curve(Y0[:,[False,False,True,False,False]], P0[:,[False,False,True,False,False]])
fpr_keras0_4, tpr_keras0_4, thresholds_kerast = roc_curve(Y0[:,[False,False,False,True,False]], P0[:,[False,False,False,True,False]])
fpr_keras0_5, tpr_keras0_5, thresholds_kerast = roc_curve(Y0[:,[False,False,False,False,True]], P0[:,[False,False,False,False,True]])
fpr_keras1_1, tpr_keras1_1, thresholds_kerast = roc_curve(Y1[:,[True,False,False,False,False]], P1[:,[True,False,False,False,False]])
fpr_keras1_2, tpr_keras1_2, thresholds_kerast = roc_curve(Y1[:,[False,True,False,False,False]], P1[:,[False,True,False,False,False]])
fpr_keras1_3, tpr_keras1_3, thresholds_kerast = roc_curve(Y1[:,[False,False,True,False,False]], P1[:,[False,False,True,False,False]])
fpr_keras1_4, tpr_keras1_4, thresholds_kerast = roc_curve(Y1[:,[False,False,False,True,False]], P1[:,[False,False,False,True,False]])
fpr_keras1_5, tpr_keras1_5, thresholds_kerast = roc_curve(Y1[:,[False,False,False,False,True]], P1[:,[False,False,False,False,True]])
np.save('EndoFull/rf_tpr_kears0_1_neg.npy', tpr_keras0_1)
np.save('EndoFull/rf_tpr_kears0_2_neg.npy', tpr_keras0_2)
np.save('EndoFull/rf_tpr_kears0_3_neg.npy', tpr_keras0_3)
np.save('EndoFull/rf_tpr_kears0_4_neg.npy', tpr_keras0_4)
np.save('EndoFull/rf_tpr_kears0_5_neg.npy', tpr_keras0_5)
np.save('EndoFull/rf_tpr_kears1_1_neg.npy', tpr_keras1_1)
np.save('EndoFull/rf_tpr_kears1_2_neg.npy', tpr_keras1_2)
np.save('EndoFull/rf_tpr_kears1_3_neg.npy', tpr_keras1_3)
np.save('EndoFull/rf_tpr_kears1_4_neg.npy', tpr_keras1_4)
np.save('EndoFull/rf_tpr_kears1_5_neg.npy', tpr_keras1_5)
np.save('EndoFull/rf_fpr_kears0_1_neg.npy', fpr_keras0_1)
np.save('EndoFull/rf_fpr_kears0_2_neg.npy', fpr_keras0_2)
np.save('EndoFull/rf_fpr_kears0_3_neg.npy', fpr_keras0_3)
np.save('EndoFull/rf_fpr_kears0_4_neg.npy', fpr_keras0_4)
np.save('EndoFull/rf_fpr_kears0_5_neg.npy', fpr_keras0_5)
np.save('EndoFull/rf_fpr_kears1_1_neg.npy', fpr_keras1_1)
np.save('EndoFull/rf_fpr_kears1_2_neg.npy', fpr_keras1_2)
np.save('EndoFull/rf_fpr_kears1_3_neg.npy', fpr_keras1_3)
np.save('EndoFull/rf_fpr_kears1_4_neg.npy', fpr_keras1_4)
np.save('EndoFull/rf_fpr_kears1_5_neg.npy', fpr_keras1_5)
pr_keras0_1, re_keras0_1, thresholds_kerast = precision_recall_curve(Y0[:,[True,False,False,False,False]], P0[:,[True,False,False,False,False]])
pr_keras0_2, re_keras0_2, thresholds_kerast = precision_recall_curve(Y0[:,[False,True,False,False,False]], P0[:,[False,True,False,False,False]])
pr_keras0_3, re_keras0_3, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,True,False,False]], P0[:,[False,False,True,False,False]])
pr_keras0_4, re_keras0_4, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,False,True,False]], P0[:,[False,False,False,True,False]])
pr_keras0_5, re_keras0_5, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,False,False,True]], P0[:,[False,False,False,False,True]])
pr_keras1_1, re_keras1_1, thresholds_kerast = precision_recall_curve(Y1[:,[True,False,False,False,False]], P1[:,[True,False,False,False,False]])
pr_keras1_2, re_keras1_2, thresholds_kerast = precision_recall_curve(Y1[:,[False,True,False,False,False]], P1[:,[False,True,False,False,False]])
pr_keras1_3, re_keras1_3, thresholds_kerast = precision_recall_curve(Y1[:,[False,False,True,False,False]], P1[:,[False,False,True,False,False]])
pr_keras1_4, re_keras1_4, thresholds_kerast = precision_recall_curve(Y1[:,[False,False,False,True,False]], P1[:,[False,False,False,True,False]])
pr_keras1_5, re_keras1_5, thresholds_kerast = precision_recall_curve(Y1[:,[False,False,False,False,True]], P1[:,[False,False,False,False,True]])
np.save('EndoFull/rf_pr_kears0_1_neg.npy', pr_keras0_1)
np.save('EndoFull/rf_pr_kears0_2_neg.npy', pr_keras0_2)
np.save('EndoFull/rf_pr_kears0_3_neg.npy', pr_keras0_3)
np.save('EndoFull/rf_pr_kears0_4_neg.npy', pr_keras0_4)
np.save('EndoFull/rf_pr_kears0_5_neg.npy', pr_keras0_5)
np.save('EndoFull/rf_pr_kears1_1_neg.npy', pr_keras1_1)
np.save('EndoFull/rf_pr_kears1_2_neg.npy', pr_keras1_2)
np.save('EndoFull/rf_pr_kears1_3_neg.npy', pr_keras1_3)
np.save('EndoFull/rf_pr_kears1_4_neg.npy', pr_keras1_4)
np.save('EndoFull/rf_pr_kears1_5_neg.npy', pr_keras1_5)
np.save('EndoFull/rf_re_kears0_1_neg.npy', re_keras0_1)
np.save('EndoFull/rf_re_kears0_2_neg.npy', re_keras0_2)
np.save('EndoFull/rf_re_kears0_3_neg.npy', re_keras0_3)
np.save('EndoFull/rf_re_kears0_4_neg.npy', re_keras0_4)
np.save('EndoFull/rf_re_kears0_5_neg.npy', re_keras0_5)
np.save('EndoFull/rf_re_kears1_1_neg.npy', re_keras1_1)
np.save('EndoFull/rf_re_kears1_2_neg.npy', re_keras1_2)
np.save('EndoFull/rf_re_kears1_3_neg.npy', re_keras1_3)
np.save('EndoFull/rf_re_kears1_4_neg.npy', re_keras1_4)
np.save('EndoFull/rf_re_kears1_5_neg.npy', re_keras1_5)

#NOW GENERATE SOME EXAMPLE PLOTS
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
plt.switch_backend('agg')

#classpred = model.predict([Xchr,Xexp], batch_size=1000)
classpred = model1.predict([Xchr], batch_size=1000)

OnP = ((classpred[:,0]>np.max(classpred[:,1:5],axis=1)))
#OnP = ((classpred:,0]>np.max(classpred[:,1:5],axis=1)))
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

np.save('EndoChr/OnTP.npy', OnTP)
np.save('EndoChr/OnFP.npy', OnFP)
np.save('EndoChr/OnTN.npy', OnTN)
np.save('EndoChr/OnFN.npy', OnFN)
np.save('EndoChr/OffTP.npy', OffTP)
np.save('EndoChr/OffFP.npy', OffFP)
np.save('EndoChr/OffTN.npy', OffTN)
np.save('EndoChr/OffFN.npy', OffFN)
np.save('EndoChr/RUTP.npy', RUTP)
np.save('EndoChr/RUFP.npy', RUFP)
np.save('EndoChr/RUTN.npy', RUTN)
np.save('EndoChr/RUFN.npy', RUFN)
np.save('EndoChr/RDTP.npy', RDTP)
np.save('EndoChr/RDFP.npy', RDFP)
np.save('EndoChr/RDTN.npy', RDTN)
np.save('EndoChr/RDFN.npy', RDFN)

xs1_1 = Xchr[OnTP,:,:]
xs2_1 = Xexp[OnTP,:,:]
ys_1 = Yalt[OnTP,:]

xs1_2 = Xchr[OffTP,:,:]
xs2_2 = Xexp[OffTP,:,:]
ys_2 = Yalt[OffTP,:]

xs1_3 = Xchr[RUTP,:,:]
xs2_3 = Xexp[RUTP,:,:]
ys_3 = Yalt[RUTP,:]

xs1_4 = Xchr[RDTP,:,:]
xs2_4 = Xexp[RDTP,:,:]
ys_4 = Yalt[RDTP,:]

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model1.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model1.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions4 = de.explain('grad*input', target_tensor*np.array([1., 0., 0., 0., 0.]), input_tensors, [xs1_1])
    attributions4B = de.explain('grad*input', target_tensor*np.array([0., 1., 0., 0., 0.]), input_tensors, [xs1_1])

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model1.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model1.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions5 = de.explain('grad*input', target_tensor*np.array([0, 0., 1., 0., 0.]), input_tensors, [xs1_2])
    attributions5B = de.explain('grad*input', target_tensor*np.array([0, 0., 0., 1., 0.]), input_tensors, [xs1_2])

np.save('EndoChr/On_attributions_OnTP.npy', attributions4[0])
np.save('EndoChr/Down_attributions_OnTP.npy', attributions4B[0])
np.save('EndoChr/Off_attributions_OffTP.npy', attributions5[0])
np.save('EndoChr/Up_attributions_OffTP.npy', attributions5B[0])

from sklearn.cluster import KMeans

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model1.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model1.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions4 = de.explain('grad*input', target_tensor*np.array([0, 0., 0., 1., 0.]), input_tensors, [xs1_3])
    attributions4B = de.explain('grad*input', target_tensor*np.array([0, 0., 1., 0., 0.]), input_tensors, [xs1_3])

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model1.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model1.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions5 = de.explain('grad*input', target_tensor*np.array([0, 1., 0., 0., 0.]), input_tensors, [xs1_4])
    attributions5B = de.explain('grad*input', target_tensor*np.array([1, 0., 0., 0., 0.]), input_tensors, [xs1_4])


np.save('EndoChr/Up_attributions_RUTP.npy', attributions4[0])
np.save('EndoChr/Off_attributions_RUTP.npy', attributions4B[0])
np.save('EndoChr/Down_attributions_RDTP.npy', attributions5[0])
np.save('EndoChr/On_attributions_RDTP.npy', attributions5B[0])

labs=["Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation"]



OnTP = np.random.permutation(OnTP)
OffTP = np.random.permutation(OffTP)
RUTP = np.random.permutation(RUTP)
RDTP = np.random.permutation(RDTP)

np.save('EndoChr/OnTP_shuffle.npy', OnTP)
np.save('EndoChr/OffTP_shuffle.npy', OffTP)
np.save('EndoChr/RUTP_shuffle.npy', RUTP)
np.save('EndoChr/RDTP_shuffle.npy', RDTP)


xs1_1 = Xchr[OnTP,:,:]
xs2_1 = Xexp[OnTP,:,:]
ys_1 = Yalt[OnTP,:]

xs1_2 = Xchr[OffTP,:,:]
xs2_2 = Xexp[OffTP,:,:]
ys_2 = Yalt[OffTP,:]

xs1_3 = Xchr[RUTP,:,:]
xs2_3 = Xexp[RUTP,:,:]
ys_3 = Yalt[RUTP,:]

xs1_4 = Xchr[RDTP,:,:]
xs2_4 = Xexp[RDTP,:,:]
ys_4 = Yalt[RDTP,:]

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions4 = de.explain('grad*input', target_tensor*np.array([1., 0., 0., 0., 0.]), input_tensors, [xs1_1,xs2_1])
    attributions4B = de.explain('grad*input', target_tensor*np.array([0., 1., 0., 0., 0.]), input_tensors, [xs1_1,xs2_1])

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions5 = de.explain('grad*input', target_tensor*np.array([0, 0., 1., 0., 0.]), input_tensors, [xs1_2,xs2_2])
    attributions5B = de.explain('grad*input', target_tensor*np.array([0, 0., 0., 1., 0.]), input_tensors, [xs1_2,xs2_2])

np.save('EndoChr/On_attributions_OnTP_perm.npy', attributions4[0])
np.save('EndoChr/Down_attributions_OnTP_perm.npy', attributions4B[0])
np.save('EndoChr/Off_attributions_OffTP_perm.npy', attributions5[0])
np.save('EndoChr/Up_attributions_OffTP_perm.npy', attributions5B[0])


with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions4 = de.explain('grad*input', target_tensor*np.array([0, 0., 0., 1., 0.]), input_tensors, [xs1_3,xs2_3])
    attributions4B = de.explain('grad*input', target_tensor*np.array([0, 0., 1., 0., 0.]), input_tensors, [xs1_3,xs2_3])

with DeepExplain(session=K.get_session()) as de:
    input_tensors = model.inputs #layers[0].input
    fModel = Model(inputs=input_tensors, outputs = model.layers[-2].output)
    target_tensor = fModel(input_tensors)
    attributions5 = de.explain('grad*input', target_tensor*np.array([0, 1., 0., 0., 0.]), input_tensors, [xs1_4,xs2_4])
    attributions5B = de.explain('grad*input', target_tensor*np.array([1, 0., 0., 0., 0.]), input_tensors, [xs1_4,xs2_4])

np.save('EndoChr/Up_attributions_RUTP_perm.npy', attributions4[0])
np.save('EndoChr/Off_attributions_RUTP_perm.npy', attributions4B[0])
np.save('EndoChr/Down_attributions_RDTP_perm.npy', attributions5[0])
np.save('EndoChr/On_attributions_RDTP_perm.npy', attributions5B[0])



A =np.squeeze(np.reshape(Xchr[trainset,:,:],(np.shape(Xchr[trainset,:,:])[0],   np.shape(Xchr[trainset,:,:])[1] *np.shape(Xchr[trainset,:,:])[2])  ))
C = A #np.concatenate((A,B),axis=1)
clf = RandomForestClassifier(n_estimators=200, max_depth=2,random_state=0)
clf.fit(C , np.argmax(Yalt[trainset,:], axis=1) )
A1 =np.squeeze(np.reshape(Xchr[valset,:,:],(np.shape(Xchr[valset,:,:])[0],   np.shape(Xchr[valset,:,:])[1] *np.shape(Xchr[testset,:,:])[2])  ))
B1 =np.squeeze(Xexp[valset,:,:] )
C1 = A1 #np.concatenate((A1,B1),axis=1)
rfpred = clf.predict_proba(C1)
np.save('EndoChr/rfpreds1.npy',rfpred)
#np.save('EndoChr/rfpredTP.npy',Yalt[testset,:])
fpr_keras, tpr_keras, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), rfpred.ravel())
pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), rfpred.ravel())
np.save('EndoChr/rf_fpr_kearss1.npy', fpr_keras)
np.save('EndoChr/rf_tpr_kearss1.npy', tpr_keras)
np.save('EndoChr/rf_pr_kearss1.npy', pr_keras)
np.save('EndoChr/rf_re_kearss1.npy', re_keras)

P0 = rfpred

clf = RandomForestClassifier(n_estimators=200, max_depth=2,random_state=1)
clf.fit(C , np.argmax(Yalt[trainset,:], axis=1) )
A1 =np.squeeze(np.reshape(Xchr[valset,:,:],(np.shape(Xchr[valset,:,:])[0],   np.shape(Xchr[valset,:,:])[1] *np.shape(Xchr[testset,:,:])[2])  ))
B1 =np.squeeze(Xexp[valset,:,:] )
C1 = A1 #np.concatenate((A1,B1),axis=1)
rfpred = clf.predict_proba(C1)
np.save('EndoChr/rfpreds2.npy',rfpred)
#np.save('EndoChr/rfpredTP.npy',Yalt[testset,:])
fpr_keras, tpr_keras, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), rfpred.ravel())
pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), rfpred.ravel())
np.save('EndoChr/rf_fpr_kearss2.npy', fpr_keras)
np.save('EndoChr/rf_tpr_kearss2.npy', tpr_keras)
np.save('EndoChr/rf_pr_kearss2.npy', pr_keras)
np.save('EndoChr/rf_re_kearss2.npy', re_keras)

P1 = rfpred
#np.save('EndoChrSomEndo/tfpr_rf.npy',fpr_keras_rf)
#np.save('EndoChrSomEndo/pr_rf.npy',pr_keras_rf)
#np.save('EndoChrSomEndo/re_rf.npy',re_keras_rf)
P1 = rfpred
Y0 = Yalt[valset,:]
Y1 = Yalt[valset,:]
fpr_keras0_1, tpr_keras0_1, thresholds_kerast = roc_curve(Y0[:,[True,False,False,False,False]], P0[:,[True,False,False,False,False]])
fpr_keras0_2, tpr_keras0_2, thresholds_kerast = roc_curve(Y0[:,[False,True,False,False,False]], P0[:,[False,True,False,False,False]])
fpr_keras0_3, tpr_keras0_3, thresholds_kerast = roc_curve(Y0[:,[False,False,True,False,False]], P0[:,[False,False,True,False,False]])
fpr_keras0_4, tpr_keras0_4, thresholds_kerast = roc_curve(Y0[:,[False,False,False,True,False]], P0[:,[False,False,False,True,False]])
fpr_keras0_5, tpr_keras0_5, thresholds_kerast = roc_curve(Y0[:,[False,False,False,False,True]], P0[:,[False,False,False,False,True]])
fpr_keras1_1, tpr_keras1_1, thresholds_kerast = roc_curve(Y1[:,[True,False,False,False,False]], P1[:,[True,False,False,False,False]])
fpr_keras1_2, tpr_keras1_2, thresholds_kerast = roc_curve(Y1[:,[False,True,False,False,False]], P1[:,[False,True,False,False,False]])
fpr_keras1_3, tpr_keras1_3, thresholds_kerast = roc_curve(Y1[:,[False,False,True,False,False]], P1[:,[False,False,True,False,False]])
fpr_keras1_4, tpr_keras1_4, thresholds_kerast = roc_curve(Y1[:,[False,False,False,True,False]], P1[:,[False,False,False,True,False]])
fpr_keras1_5, tpr_keras1_5, thresholds_kerast = roc_curve(Y1[:,[False,False,False,False,True]], P1[:,[False,False,False,False,True]])
np.save('EndoChr/rf_tpr_kears0_1.npy', tpr_keras0_1)
np.save('EndoChr/rf_tpr_kears0_2.npy', tpr_keras0_2)
np.save('EndoChr/rf_tpr_kears0_3.npy', tpr_keras0_3)
np.save('EndoChr/rf_tpr_kears0_4.npy', tpr_keras0_4)
np.save('EndoChr/rf_tpr_kears0_5.npy', tpr_keras0_5)
np.save('EndoChr/rf_tpr_kears1_1.npy', tpr_keras1_1)
np.save('EndoChr/rf_tpr_kears1_2.npy', tpr_keras1_2)
np.save('EndoChr/rf_tpr_kears1_3.npy', tpr_keras1_3)
np.save('EndoChr/rf_tpr_kears1_4.npy', tpr_keras1_4)
np.save('EndoChr/rf_tpr_kears1_5.npy', tpr_keras1_5)
np.save('EndoChr/rf_fpr_kears0_1.npy', fpr_keras0_1)
np.save('EndoChr/rf_fpr_kears0_2.npy', fpr_keras0_2)
np.save('EndoChr/rf_fpr_kears0_3.npy', fpr_keras0_3)
np.save('EndoChr/rf_fpr_kears0_4.npy', fpr_keras0_4)
np.save('EndoChr/rf_fpr_kears0_5.npy', fpr_keras0_5)
np.save('EndoChr/rf_fpr_kears1_1.npy', fpr_keras1_1)
np.save('EndoChr/rf_fpr_kears1_2.npy', fpr_keras1_2)
np.save('EndoChr/rf_fpr_kears1_3.npy', fpr_keras1_3)
np.save('EndoChr/rf_fpr_kears1_4.npy', fpr_keras1_4)
np.save('EndoChr/rf_fpr_kears1_5.npy', fpr_keras1_5)
pr_keras0_1, re_keras0_1, thresholds_kerast = precision_recall_curve(Y0[:,[True,False,False,False,False]], P0[:,[True,False,False,False,False]])
pr_keras0_2, re_keras0_2, thresholds_kerast = precision_recall_curve(Y0[:,[False,True,False,False,False]], P0[:,[False,True,False,False,False]])
pr_keras0_3, re_keras0_3, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,True,False,False]], P0[:,[False,False,True,False,False]])
pr_keras0_4, re_keras0_4, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,False,True,False]], P0[:,[False,False,False,True,False]])
pr_keras0_5, re_keras0_5, thresholds_kerast = precision_recall_curve(Y0[:,[False,False,False,False,True]], P0[:,[False,False,False,False,True]])
pr_keras1_1, re_keras1_1, thresholds_kerast = precision_recall_curve(Y1[:,[True,False,False,False,False]], P1[:,[True,False,False,False,False]])
pr_keras1_2, re_keras1_2, thresholds_kerast = precision_recall_curve(Y1[:,[False,True,False,False,False]], P1[:,[False,True,False,False,False]])
pr_keras1_3, re_keras1_3, thresholds_kerast = precision_recall_curve(Y1[:,[False,False,True,False,False]], P1[:,[False,False,True,False,False]])
pr_keras1_4, re_keras1_4, thresholds_kerast = precision_recall_curve(Y1[:,[False,False,False,True,False]], P1[:,[False,False,False,True,False]])
pr_keras1_5, re_keras1_5, thresholds_kerast = precision_recall_curve(Y1[:,[False,False,False,False,True]], P1[:,[False,False,False,False,True]])
np.save('EndoChr/rf_pr_kears0_1.npy', pr_keras0_1)
np.save('EndoChr/rf_pr_kears0_2.npy', pr_keras0_2)
np.save('EndoChr/rf_pr_kears0_3.npy', pr_keras0_3)
np.save('EndoChr/rf_pr_kears0_4.npy', pr_keras0_4)
np.save('EndoChr/rf_pr_kears0_5.npy', pr_keras0_5)
np.save('EndoChr/rf_pr_kears1_1.npy', pr_keras1_1)
np.save('EndoChr/rf_pr_kears1_2.npy', pr_keras1_2)
np.save('EndoChr/rf_pr_kears1_3.npy', pr_keras1_3)
np.save('EndoChr/rf_pr_kears1_4.npy', pr_keras1_4)
np.save('EndoChr/rf_pr_kears1_5.npy', pr_keras1_5)
np.save('EndoChr/rf_re_kears0_1.npy', re_keras0_1)
np.save('EndoChr/rf_re_kears0_2.npy', re_keras0_2)
np.save('EndoChr/rf_re_kears0_3.npy', re_keras0_3)
np.save('EndoChr/rf_re_kears0_4.npy', re_keras0_4)
np.save('EndoChr/rf_re_kears0_5.npy', re_keras0_5)
np.save('EndoChr/rf_re_kears1_1.npy', re_keras1_1)
np.save('EndoChr/rf_re_kears1_2.npy', re_keras1_2)
np.save('EndoChr/rf_re_kears1_3.npy', re_keras1_3)
np.save('EndoChr/rf_re_kears1_4.npy', re_keras1_4)
np.save('EndoChr/rf_re_kears1_5.npy', re_keras1_5)


#from sklearn.ensemble import RandomForestClassifier
#clf = RandomForestClassifier(n_estimators=100, max_depth=2,random_state=0)
#clf.fit([np.squeeze(np.reshape(Xchr[trainset,:,:],(np.shape(Xchr[trainset,:,:])[0],   np.shape(Xchr[trainset,:,:])[1] *np.shape(Xchr[trainset,:,:])[2])  )) , np.argmax(Yalt[trainset,:], axis=1),Xexp[trainset,:,:]))
#rfpred = clf.predict_proba(np.squeeze(np.reshape(Xchr[testset,:,:],(np.shape(Xchr[testset,:,:])[0],   np.shape(Xchr[testset,:,:])[1] *np.shape(Xchr[testset,:,:])[2])  )))
#fpr_keras_rf, tpr_keras_rf, thresholds_keras = roc_curve(Yalt[testset,:].ravel(), rfpred.ravel())
#pr_keras_rf, re_keras_rf, thresholds_keras = precision_recall_curve(Yalt[testset,:].ravel(), rfpred.ravel())
#np.save('EndoChrSomEndo/tpr_rf.npy',tpr_keras_rf)
#np.save('EndoChrSomEndo/tfpr_rf.npy',fpr_keras_rf)
#np.save('EndoChrSomEndo/pr_rf.npy',pr_keras_rf)
#np.save('EndoChrSomEndo/re_rf.npy',re_keras_rf)
