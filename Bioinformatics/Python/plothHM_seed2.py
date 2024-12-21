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
model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model1 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_neg_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model2 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_pos_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)

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

np.save('EndoFull/OnTP_seed=2.npy', OnTP)
np.save('EndoFull/OnFP_seed=2.npy', OnFP)
np.save('EndoFull/OnTN_seed=2.npy', OnTN)
np.save('EndoFull/OnFN_seed=2.npy', OnFN)
np.save('EndoFull/OffTP_seed=2.npy', OffTP)
np.save('EndoFull/OffFP_seed=2.npy', OffFP)
np.save('EndoFull/OffTN_seed=2.npy', OffTN)
np.save('EndoFull/OffFN_seed=2.npy', OffFN)
np.save('EndoFull/RUTP_seed=2.npy', RUTP)
np.save('EndoFull/RUFP_seed=2.npy', RUFP)
np.save('EndoFull/RUTN_seed=2.npy', RUTN)
np.save('EndoFull/RUFN_seed=2.npy', RUFN)
np.save('EndoFull/RDTP_seed=2.npy', RDTP)
np.save('EndoFull/RDFP_seed=2.npy', RDFP)
np.save('EndoFull/RDTN_seed=2.npy', RDTN)
np.save('EndoFull/RDFN_seed=2.npy', RDFN)


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

np.save('EndoFull/Up_attributions_RUTP_seed=2.npy', attributions4)
np.save('EndoFull/Off_attributions_RUTP_seed=2.npy', attributions4B)
np.save('EndoFull/Down_attributions_RDTP_seed=2.npy', attributions5)
np.save('EndoFull/On_attributions_RDTP_seed=2.npy', attributions5B)


#from sklearn.cluster import KMeans

kmeans1 = KMeans(n_clusters=5, random_state=0).fit(np.reshape(attributions4[0],(np.shape(attributions4[0])[0],np.shape(attributions4[0])[1]*np.shape(attributions4[0])[2])))
kmeans1.labels_

hm1 = np.mean(attributions4[0],axis=0)
hm2 = np.mean(attributions4[0][kmeans1.labels_==0,:,:],axis=0)
hm3 = np.mean(attributions4[0][kmeans1.labels_==1,:,:],axis=0)
hm4 = np.mean(attributions4[0][kmeans1.labels_==2,:,:],axis=0)
hm5 = np.mean(attributions4[0][kmeans1.labels_==3,:,:],axis=0)
hm6 = np.mean(attributions4[0][kmeans1.labels_==4,:,:],axis=0)

kmeans2 = KMeans(n_clusters=5, random_state=0).fit(np.reshape(attributions5[0],(np.shape(attributions5[0])[0],np.shape(attributions5[0])[1]*np.shape(attributions5[0])[2])))
kmeans2.labels_

hm7 = np.mean(attributions5[0],axis=0)
hm8 = np.mean(attributions5[0][kmeans2.labels_==0,:,:],axis=0)
hm9 = np.mean(attributions5[0][kmeans2.labels_==1,:,:],axis=0)
hm10 = np.mean(attributions5[0][kmeans2.labels_==2,:,:],axis=0)
hm11 = np.mean(attributions5[0][kmeans2.labels_==3,:,:],axis=0)
hm12 = np.mean(attributions5[0][kmeans2.labels_==4,:,:],axis=0)


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


np.save('EndoFull/Up_attributions_RUTP_seed=2.npy', attributions4)
np.save('EndoFull/Off_attributions_RUTP_seed=2.npy', attributions4B)
np.save('EndoFull/Down_attributions_RDTP_seed=2.npy', attributions5)
np.save('EndoFull/On_attributions_RDTP_seed=2.npy', attributions5B)


kmeans3 = KMeans(n_clusters=5, random_state=0).fit(np.reshape(attributions4[0],(np.shape(attributions4[0])[0],np.shape(attributions4[0])[1]*np.shape(attributions4[0])[2])))
kmeans3.labels_

hm13 = np.mean(attributions4[0],axis=0)
hm14 = np.mean(attributions4[0][kmeans3.labels_==0,:,:],axis=0)
hm15 = np.mean(attributions4[0][kmeans3.labels_==1,:,:],axis=0)
hm16 = np.mean(attributions4[0][kmeans3.labels_==2,:,:],axis=0)
hm17 = np.mean(attributions4[0][kmeans3.labels_==3,:,:],axis=0)
hm18 = np.mean(attributions4[0][kmeans3.labels_==4,:,:],axis=0)


kmeans4 = KMeans(n_clusters=5, random_state=0).fit(np.reshape(attributions5[0],(np.shape(attributions5[0])[0],np.shape(attributions5[0])[1]*np.shape(attributions5[0])[2])))
kmeans4.labels_

hm19 = np.mean(attributions5[0],axis=0)
hm20 = np.mean(attributions5[0][kmeans4.labels_==0,:,:],axis=0)
hm21 = np.mean(attributions5[0][kmeans4.labels_==1,:,:],axis=0)
hm22 = np.mean(attributions5[0][kmeans4.labels_==2,:,:],axis=0)
hm23 = np.mean(attributions5[0][kmeans4.labels_==3,:,:],axis=0)
hm24 = np.mean(attributions5[0][kmeans4.labels_==4,:,:],axis=0)




np.save('OnTP_ac_deeper_seed=2.npy', hm1)
np.save('OnTP_ac_deeper_cl1_seed=2.npy', hm2)
np.save('OnTP_ac_deeper_cl2_seed=2.npy', hm3)
np.save('OnTP_ac_deeper_cl3_seed=2.npy', hm4)
np.save('OnTP_ac_deeper_cl4_seed=2.npy', hm5)
np.save('OnTP_ac_deeper_cl5_seed=2.npy', hm6)

np.save('OffTP_ac_deeper_seed=2.npy', hm7)
np.save('OffTP_ac_deeper_cl1_seed=2.npy', hm8)
np.save('OffTP_ac_deeper_cl2_seed=2.npy', hm9)
np.save('OffTP_ac_deeper_cl3_seed=2.npy', hm10)
np.save('OffTP_ac_deeper_cl4_seed=2.npy', hm11)
np.save('OffTP_ac_deeper_cl5_seed=2.npy', hm12)

np.save('kmeans1_seed=2.npy', kmeans1)


np.savetxt('OnTP_ac_deeper_seed=2.csv', hm1, delimiter=",")
np.savetxt('OnTP_ac_deeper_cl1_seed=2.csv', hm2, delimiter=",")
np.savetxt('OnTP_ac_deeper_cl2_seed=2.csv', hm3, delimiter=",")
np.savetxt('OnTP_ac_deeper_cl3_seed=2.csv', hm4, delimiter=",")
np.savetxt('OnTP_ac_deeper_cl4_seed=2.csv', hm5, delimiter=",")
np.savetxt('OnTP_ac_deeper_cl5_seed=2.csv', hm6, delimiter=",")

np.savetxt('OffTP_ac_deeper_seed=2.csv', hm7, delimiter=",")
np.savetxt('OffTP_ac_deeper_cl1_seed=2.csv', hm8, delimiter=",")
np.savetxt('OffTP_ac_deeper_cl2_seed=2.csv', hm9, delimiter=",")
np.savetxt('OffTP_ac_deeper_cl3_seed=2.csv', hm10, delimiter=",")
np.savetxt('OffTP_ac_deeper_cl4_seed=2.csv', hm11, delimiter=",")
np.savetxt('OffTP_ac_deeper_cl5_seed=2.csv', hm12, delimiter=",")

#csv_matrix = np.savetxt("kmeans1.csv", kmeans1, delimiter=",")
#csv_matrix = np.savetxt("kmeans1.csv", kmeans1, delimiter=",")


labs=["Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation"]

#NOW GENERATE SOME EXAMPLE PLOTS
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
plt.switch_backend('agg')

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm1), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.008)
plt.colorbar()
plt.savefig('EndoFull/OnTP_ac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm2), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.008)
plt.colorbar()
plt.savefig('EndoFull/ONTP_ac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm3), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.008)
plt.colorbar()
plt.savefig('EndoFull/ONTP_ac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm4), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.008)
plt.colorbar()
plt.savefig('EndoFull/ONTP_ac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm5), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.008)
plt.colorbar()
plt.savefig('EndoFull/ONTP_ac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm6), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.008)
plt.colorbar()
plt.savefig('EndoFull/ONTP_ac_cl5_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm7), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.005)
plt.colorbar()
plt.savefig('EndoFull/OffTP_ac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm8), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.005)
plt.colorbar()
plt.savefig('EndoFull/OffTP_ac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm9), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.005)
plt.colorbar()
plt.savefig('EndoFull/OffTP_ac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm10), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.005)
plt.colorbar()
plt.savefig('EndoFull/OffTP_ac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm11), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.005)
plt.colorbar()
plt.savefig('EndoFull/OffTP_ac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm12), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)plt.clim(-0.002,0.005)
plt.clim(-0.002,0.005)
plt.colorbar()
plt.savefig('EndoFull/OffTP_ac_cl5_deeper_seed=2.pdf')
plt.close()


plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm13), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.008)

plt.colorbar()
plt.savefig('EndoFull/RUTP_ac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm14), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.colorbar()
plt.savefig('EndoFull/RUTP_ac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm15), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.008)

plt.colorbar()
plt.savefig('EndoFull/RUTP_ac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm16), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.008)

plt.colorbar()
plt.savefig('EndoFull/RUTP_ac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm17), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.008)

plt.colorbar()
plt.savefig('EndoFull/RUTP_ac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm18), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.008)

plt.colorbar()
plt.savefig('EndoFull/RUTP_ac_cl5_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm19), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFull/RDTP_ac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm20), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFull/RDTP_ac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm21), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.colorbar()
plt.savefig('EndoFull/RDTP_ac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm22), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
#plt.clim(-0.002,0.008)
plt.colorbar()
plt.savefig('EndoFull/RDTP_ac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm23), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
#plt.clim(-0.002,0.008)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFull/RDTP_ac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm24), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFull/RDTP_ac_cl5_deeper_seed=2.pdf')
plt.close()

np.save('EndoFull/hm1_delta_seed=2.npy', hm1)
np.save('EndoFull/hm2_delta_seed=2.npy', hm2)
np.save('EndoFull/hm3_delta_seed=2.npy', hm3)
np.save('EndoFull/hm4_delta_seed=2.npy', hm4)
np.save('EndoFull/hm5_delta_seed=2.npy', hm5)
np.save('EndoFull/hm6_delta_seed=2.npy', hm6)
np.save('EndoFull/hm7_delta_seed=2.npy', hm7)
np.save('EndoFull/hm8_delta_seed=2.npy', hm8)
np.save('EndoFull/hm9_delta_seed=2.npy', hm9)
np.save('EndoFull/hm10_delta_seed=2.npy', hm10)
np.save('EndoFull/hm11_delta_seed=2.npy', hm11)
np.save('EndoFull/hm12_delta_seed=2.npy', hm12)
np.save('EndoFull/hm13_delta_seed=2.npy', hm13)
np.save('EndoFull/hm14_delta_seed=2.npy', hm14)
np.save('EndoFull/hm15_delta_seed=2.npy', hm15)
np.save('EndoFull/hm16_delta_seed=2.npy', hm16)
np.save('EndoFull/hm17_delta_seed=2.npy', hm17)
np.save('EndoFull/hm18_delta_seed=2.npy', hm18)
np.save('EndoFull/hm19_delta_seed=2.npy', hm19)
np.save('EndoFull/hm20_delta_seed=2.npy', hm20)
np.save('EndoFull/hm21_delta_seed=2.npy', hm21)
np.save('EndoFull/hm22_delta_seed=2.npy', hm22)
np.save('EndoFull/hm23_delta_seed=2.npy', hm23)
np.save('EndoFull/hm24_delta_seed=2.npy', hm24)



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


kmeans1 = KMeans(n_clusters=5, random_state=0).fit( np.reshape(attributions4[0],(np.shape(attributions4[0])[0],np.shape(attributions4[0])[1]*np.shape(attributions4[0])[2]))-np.reshape(attributions4B[0],(np.shape(attributions4B[0])[0],np.shape(attributions4B[0])[1]*np.shape(attributions4B[0])[2])) )
kmeans1.labels_

hm1 = np.mean(attributions4[0]-attributions4B[0],axis=0)
hm2 = np.mean(attributions4[0][kmeans1.labels_==0,:,:]-attributions4B[0][kmeans1.labels_==0,:,:],axis=0)
hm3 = np.mean(attributions4[0][kmeans1.labels_==1,:,:]-attributions4B[0][kmeans1.labels_==1,:,:],axis=0)
hm4 = np.mean(attributions4[0][kmeans1.labels_==2,:,:]-attributions4B[0][kmeans1.labels_==2,:,:],axis=0)
hm5 = np.mean(attributions4[0][kmeans1.labels_==3,:,:]-attributions4B[0][kmeans1.labels_==3,:,:],axis=0)
hm6 = np.mean(attributions4[0][kmeans1.labels_==4,:,:]-attributions4B[0][kmeans1.labels_==4,:,:],axis=0)

kmeans2 = KMeans(n_clusters=5, random_state=0).fit(
                                                   np.reshape(attributions5[0],(np.shape(attributions5[0])[0],np.shape(attributions5[0])[1]*np.shape(attributions5[0])[2]))-np.reshape(attributions5B[0],(np.shape(attributions5B[0])[0],np.shape(attributions5B[0])[1]*np.shape(attributions5B[0])[2])) )
kmeans2.labels_

hm7 = np.mean(attributions5[0]-attributions5B[0],axis=0)
hm8 = np.mean(attributions5[0][kmeans2.labels_==0,:,:]-attributions5B[0][kmeans2.labels_==0,:,:],axis=0)
hm9 = np.mean(attributions5[0][kmeans2.labels_==1,:,:]-attributions5B[0][kmeans2.labels_==1,:,:],axis=0)
hm10 = np.mean(attributions5[0][kmeans2.labels_==2,:,:]-attributions5B[0][kmeans2.labels_==2,:,:],axis=0)
hm11 = np.mean(attributions5[0][kmeans2.labels_==3,:,:]-attributions5B[0][kmeans2.labels_==3,:,:],axis=0)
hm12 = np.mean(attributions5[0][kmeans2.labels_==4,:,:]-attributions5B[0][kmeans2.labels_==4,:,:],axis=0)


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



kmeans3 = KMeans(n_clusters=5, random_state=0).fit(
                                                   np.reshape(attributions4[0],(np.shape(attributions4[0])[0],np.shape(attributions4[0])[1]*np.shape(attributions4[0])[2]))- np.reshape(attributions4B[0],(np.shape(attributions4B[0])[0],np.shape(attributions4B[0])[1]*np.shape(attributions4B[0])[2])) )
kmeans3.labels_

hm13 = np.mean(attributions4[0]-attributions4B[0],axis=0)
hm14 = np.mean(attributions4[0][kmeans3.labels_==0,:,:]-attributions4B[0][kmeans3.labels_==0,:,:],axis=0)
hm15 = np.mean(attributions4[0][kmeans3.labels_==1,:,:]-attributions4B[0][kmeans3.labels_==1,:,:],axis=0)
hm16 = np.mean(attributions4[0][kmeans3.labels_==2,:,:]-attributions4B[0][kmeans3.labels_==2,:,:],axis=0)
hm17 = np.mean(attributions4[0][kmeans3.labels_==3,:,:]-attributions4B[0][kmeans3.labels_==3,:,:],axis=0)
hm18 = np.mean(attributions4[0][kmeans3.labels_==4,:,:]-attributions4B[0][kmeans3.labels_==4,:,:],axis=0)


kmeans4 = KMeans(n_clusters=5, random_state=0).fit(
                                                   np.reshape(attributions5[0],(np.shape(attributions5[0])[0],np.shape(attributions5[0])[1]*np.shape(attributions5[0])[2]))-np.reshape(attributions5B[0],(np.shape(attributions5B[0])[0],np.shape(attributions5B[0])[1]*np.shape(attributions5B[0])[2])) )
kmeans4.labels_

hm19 = np.mean(attributions5[0]-attributions5B[0],axis=0)
hm20 = np.mean(attributions5[0][kmeans4.labels_==0,:,:]-attributions5B[0][kmeans4.labels_==0,:,:],axis=0)
hm21 = np.mean(attributions5[0][kmeans4.labels_==1,:,:]-attributions5B[0][kmeans4.labels_==1,:,:],axis=0)
hm22 = np.mean(attributions5[0][kmeans4.labels_==2,:,:]-attributions5B[0][kmeans4.labels_==2,:,:],axis=0)
hm23 = np.mean(attributions5[0][kmeans4.labels_==3,:,:]-attributions5B[0][kmeans4.labels_==3,:,:],axis=0)
hm24 = np.mean(attributions5[0][kmeans4.labels_==4,:,:]-attributions5B[0][kmeans4.labels_==4,:,:],axis=0)

labs=["Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation"]


plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm1), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.001,0.002)
plt.colorbar()
plt.savefig('EndoFull/OnTP_Deltac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm2), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
#plt.clim(-0.005,0.002)
plt.clim(-0.001,0.002)
plt.colorbar()
plt.savefig('EndoFull/ONTP_Deltac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm3), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
#plt.clim(-0.005,0.002)
plt.clim(-0.001,0.002)

plt.colorbar()
plt.savefig('EndoFull/ONTP_Deltac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm4), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
#plt.clim(-0.005,0.002)
plt.clim(-0.001,0.002)

plt.colorbar()
plt.savefig('EndoFull/ONTP_Deltac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm5), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
#plt.clim(-0.005,0.002)
plt.clim(-0.001,0.002)

plt.colorbar()
plt.savefig('EndoFull/ONTP_Deltac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm6), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
#plt.clim(-0.005,0.002)
plt.clim(-0.001,0.002)
plt.colorbar()
plt.savefig('EndoFull/ONTP_Deltac_cl5_deeper_seed=2.pdf')
plt.close()


plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm7), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)
plt.colorbar()
plt.savefig('EndoFull/OffTP_Deltac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm8), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)

plt.colorbar()
plt.savefig('EndoFull/OffTP_Deltac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm9), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)

plt.colorbar()
plt.savefig('EndoFull/OffTP_Deltac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm10), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)

plt.colorbar()
plt.savefig('EndoFull/OffTP_Deltac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm11), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)

plt.colorbar()
plt.savefig('EndoFull/OffTP_Deltac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm12), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)
plt.colorbar()
plt.savefig('EndoFull/OffTP_Deltac_cl5_deeper_seed=2.pdf')
plt.close()






plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm13), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.01)
plt.colorbar()
plt.savefig('EndoFull/RUTP_Deltac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm14), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.01)
plt.colorbar()
plt.savefig('EndoFull/RUTP_Deltac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm15), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.01)
plt.colorbar()
plt.savefig('EndoFull/RUTP_Deltac_cl2_deeper_seed=2.pdf')
plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm16), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.01)
plt.colorbar()
plt.savefig('EndoFull/RUTP_Deltac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm17), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.01)
plt.colorbar()
plt.savefig('EndoFull/RUTP_Deltac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm18), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.01)
plt.colorbar()
plt.savefig('EndoFull/RUTP_Deltac_cl5_deeper_seed=2.pdf')
plt.close()



plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm19), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.001,0.01)
plt.colorbar()
plt.savefig('EndoFull/RDTP_Deltac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm20), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.001,0.01)
plt.colorbar()
plt.savefig('EndoFull/RDTP_Deltac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm21), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.001,0.01)
plt.colorbar()
plt.savefig('EndoFull/RDTP_Deltac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm22), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.001,0.01)
plt.colorbar()
plt.savefig('EndoFull/RDTP_Deltac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm23), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.001,0.01)
plt.colorbar()
plt.savefig('EndoFull/RDTP_Deltac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm24), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.001,0.01)
plt.colorbar()
plt.savefig('EndoFull/RDTP_Deltac_cl5_deeper_seed=2.pdf')
plt.close()


np.save('EndoFull/hm1_delta_seed=2.npy', hm1)
np.save('EndoFull/hm2_delta_seed=2.npy', hm2)
np.save('EndoFull/hm3_delta_seed=2.npy', hm3)
np.save('EndoFull/hm4_delta_seed=2.npy', hm4)
np.save('EndoFull/hm5_delta_seed=2.npy', hm5)
np.save('EndoFull/hm6_delta_seed=2.npy', hm6)
np.save('EndoFull/hm7_delta_seed=2.npy', hm7)
np.save('EndoFull/hm8_delta_seed=2.npy', hm8)
np.save('EndoFull/hm9_delta_seed=2.npy', hm9)
np.save('EndoFull/hm10_delta_seed=2.npy', hm10)
np.save('EndoFull/hm11_delta_seed=2.npy', hm11)
np.save('EndoFull/hm12_delta_seed=2.npy', hm12)
np.save('EndoFull/hm13_delta_seed=2.npy', hm13)
np.save('EndoFull/hm14_delta_seed=2.npy', hm14)
np.save('EndoFull/hm15_delta_seed=2.npy', hm15)
np.save('EndoFull/hm16_delta_seed=2.npy', hm16)
np.save('EndoFull/hm17_delta_seed=2.npy', hm17)
np.save('EndoFull/hm18_delta_seed=2.npy', hm18)
np.save('EndoFull/hm19_delta_seed=2.npy', hm19)
np.save('EndoFull/hm20_delta_seed=2.npy', hm20)
np.save('EndoFull/hm21_delta_seed=2.npy', hm21)
np.save('EndoFull/hm22_delta_seed=2.npy', hm22)
np.save('EndoFull/hm23_delta_seed=2.npy', hm23)
np.save('EndoFull/hm24_delta_seed=2.npy', hm24)




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

np.save('EndoChr/OnTP_seed=2.npy', OnTP)
np.save('EndoChr/OnFP_seed=2.npy', OnFP)
np.save('EndoChr/OnTN_seed=2.npy', OnTN)
np.save('EndoChr/OnFN_seed=2.npy', OnFN)
np.save('EndoChr/OffTP_seed=2.npy', OffTP)
np.save('EndoChr/OffFP_seed=2.npy', OffFP)
np.save('EndoChr/OffTN_seed=2.npy', OffTN)
np.save('EndoChr/OffFN_seed=2.npy', OffFN)
np.save('EndoChr/RUTP_seed=2.npy', RUTP)
np.save('EndoChr/RUFP_seed=2.npy', RUFP)
np.save('EndoChr/RUTN_seed=2.npy', RUTN)
np.save('EndoChr/RUFN_seed=2.npy', RUFN)
np.save('EndoChr/RDTP_seed=2.npy', RDTP)
np.save('EndoChr/RDFP_seed=2.npy', RDFP)
np.save('EndoChr/RDTN_seed=2.npy', RDTN)
np.save('EndoChr/RDFN_seed=2.npy', RDFN)

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


np.save('EndoChr/On_attributions_OnTP_seed=2.npy', attributions4)
np.save('EndoChr/Down_attributions_OnTP_seed=2.npy', attributions4B)
np.save('EndoChr/Off_attributions_OffTP_seed=2.npy', attributions5)
np.save('EndoChr/Up_attributions_OffTP_seed=2.npy', attributions5B)


from sklearn.cluster import KMeans

kmeans1 = KMeans(n_clusters=5, random_state=0).fit(np.reshape(attributions4[0],(np.shape(attributions4[0])[0],np.shape(attributions4[0])[1]*np.shape(attributions4[0])[2])))
kmeans1.labels_

hm1 = np.mean(attributions4[0],axis=0)
hm2 = np.mean(attributions4[0][kmeans1.labels_==0,:,:],axis=0)
hm3 = np.mean(attributions4[0][kmeans1.labels_==1,:,:],axis=0)
hm4 = np.mean(attributions4[0][kmeans1.labels_==2,:,:],axis=0)
hm5 = np.mean(attributions4[0][kmeans1.labels_==3,:,:],axis=0)
hm6 = np.mean(attributions4[0][kmeans1.labels_==4,:,:],axis=0)

kmeans2 = KMeans(n_clusters=5, random_state=0).fit(np.reshape(attributions5[0],(np.shape(attributions5[0])[0],np.shape(attributions5[0])[1]*np.shape(attributions5[0])[2])))
kmeans2.labels_

hm7 = np.mean(attributions5[0],axis=0)
hm8 = np.mean(attributions5[0][kmeans2.labels_==0,:,:],axis=0)
hm9 = np.mean(attributions5[0][kmeans2.labels_==1,:,:],axis=0)
hm10 = np.mean(attributions5[0][kmeans2.labels_==2,:,:],axis=0)
hm11 = np.mean(attributions5[0][kmeans2.labels_==3,:,:],axis=0)
hm12 = np.mean(attributions5[0][kmeans2.labels_==4,:,:],axis=0)


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

np.save('EndoChr/Up_attributions_RUTP_seed=2.npy', attributions4)
np.save('EndoChr/Off_attributions_RUTP_seed=2.npy', attributions4B)
np.save('EndoChr/Down_attributions_RDTP_seed=2.npy', attributions5)
np.save('EndoChr/On_attributions_RDTP_seed=2.npy', attributions5B)


kmeans3 = KMeans(n_clusters=5, random_state=0).fit(np.reshape(attributions4[0],(np.shape(attributions4[0])[0],np.shape(attributions4[0])[1]*np.shape(attributions4[0])[2])))
kmeans3.labels_

hm13 = np.mean(attributions4[0],axis=0)
hm14 = np.mean(attributions4[0][kmeans3.labels_==0,:,:],axis=0)
hm15 = np.mean(attributions4[0][kmeans3.labels_==1,:,:],axis=0)
hm16 = np.mean(attributions4[0][kmeans3.labels_==2,:,:],axis=0)
hm17 = np.mean(attributions4[0][kmeans3.labels_==3,:,:],axis=0)
hm18 = np.mean(attributions4[0][kmeans3.labels_==4,:,:],axis=0)


kmeans4 = KMeans(n_clusters=5, random_state=0).fit(np.reshape(attributions5[0],(np.shape(attributions5[0])[0],np.shape(attributions5[0])[1]*np.shape(attributions5[0])[2])))
kmeans4.labels_

hm19 = np.mean(attributions5[0],axis=0)
hm20 = np.mean(attributions5[0][kmeans4.labels_==0,:,:],axis=0)
hm21 = np.mean(attributions5[0][kmeans4.labels_==1,:,:],axis=0)
hm22 = np.mean(attributions5[0][kmeans4.labels_==2,:,:],axis=0)
hm23 = np.mean(attributions5[0][kmeans4.labels_==3,:,:],axis=0)
hm24 = np.mean(attributions5[0][kmeans4.labels_==4,:,:],axis=0)


labs=["Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation"]

#NOW GENERATE SOME EXAMPLE PLOTS
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
plt.switch_backend('agg')

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm1), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.001,0.005)
plt.colorbar()
plt.savefig('EndoChr/OnTP_ac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm2), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.005)
plt.colorbar()
plt.savefig('EndoChr/ONTP_ac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm3), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.005)
plt.colorbar()
plt.savefig('EndoChr/ONTP_ac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm4), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.005)
plt.colorbar()
plt.savefig('EndoChr/ONTP_ac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm5), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.005)
plt.colorbar()
plt.savefig('EndoChr/ONTP_ac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm6), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.005)
plt.colorbar()
plt.savefig('EndoChr/ONTP_ac_cl5_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm7), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_ac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm8), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_ac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm9), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_ac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm10), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_ac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm11), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_ac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm12), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_ac_cl5_deeper_seed=2.pdf')
plt.close()




plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm13), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.0,0.004)
plt.colorbar()
plt.savefig('EndoChr/RUTP_ac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm14), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0,0.004)
plt.colorbar()
plt.savefig('EndoChr/RUTP_ac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm15), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0,0.004)
plt.colorbar()
plt.savefig('EndoChr/RUTP_ac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm16), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0,0.004)
plt.colorbar()
plt.savefig('EndoChr/RUTP_ac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm17), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0,0.004)
plt.colorbar()
plt.savefig('EndoChr/RUTP_ac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm18), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0,0.004)
plt.colorbar()
plt.savefig('EndoChr/RUTP_ac_cl5_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm19), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.0006,0.001)
plt.colorbar()
plt.savefig('EndoChr/RDTP_ac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm20), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0006,0.001)
plt.colorbar()
plt.savefig('EndoChr/RDTP_ac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm21), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0006,0.001)
plt.colorbar()
plt.savefig('EndoChr/RDTP_ac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm22), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0006,0.001)
plt.colorbar()
plt.savefig('EndoChr/RDTP_ac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm23), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.colorbar()
plt.savefig('EndoChr/RDTP_ac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm24), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0006,0.001)
plt.colorbar()
plt.savefig('EndoChr/RDTP_ac_cl5_deeper_seed=2.pdf')
plt.close()





np.save('EndoChr/hm1_seed=2.npy', hm1)
np.save('EndoChr/hm2_seed=2.npy', hm2)
np.save('EndoChr/hm3_seed=2.npy', hm3)
np.save('EndoChr/hm4_seed=2.npy', hm4)
np.save('EndoChr/hm5_seed=2.npy', hm4)
np.save('EndoChr/hm6_seed=2.npy', hm6)
np.save('EndoChr/hm7_seed=2.npy', hm7)
np.save('EndoChr/hm8_seed=2.npy', hm8)
np.save('EndoChr/hm9_seed=2.npy', hm9)
np.save('EndoChr/hm10_seed=2.npy', hm10)
np.save('EndoChr/hm11_seed=2.npy', hm11)
np.save('EndoChr/hm12_seed=2.npy', hm12)
np.save('EndoChr/hm13_seed=2.npy', hm13)
np.save('EndoChr/hm14_seed=2.npy', hm14)
np.save('EndoChr/hm15_seed=2.npy', hm15)
np.save('EndoChr/hm16_seed=2.npy', hm16)
np.save('EndoChr/hm17_seed=2.npy', hm17)
np.save('EndoChr/hm18_seed=2.npy', hm18)
np.save('EndoChr/hm19_seed=2.npy', hm19)
np.save('EndoChr/hm20_seed=2.npy', hm20)
np.save('EndoChr/hm21_seed=2.npy', hm21)
np.save('EndoChr/hm22_seed=2.npy', hm22)
np.save('EndoChr/hm23_seed=2.npy', hm23)
np.save('EndoChr/hm24_seed=2.npy', hm24)




#from sklearn.cluster import KMeans

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


kmeans1 = KMeans(n_clusters=5, random_state=0).fit( np.reshape(attributions4[0],(np.shape(attributions4[0])[0],np.shape(attributions4[0])[1]*np.shape(attributions4[0])[2]))-np.reshape(attributions4B[0],(np.shape(attributions4B[0])[0],np.shape(attributions4B[0])[1]*np.shape(attributions4B[0])[2])) )
kmeans1.labels_

hm1 = np.mean(attributions4[0]-attributions4B[0],axis=0)
hm2 = np.mean(attributions4[0][kmeans1.labels_==0,:,:]-attributions4B[0][kmeans1.labels_==0,:,:],axis=0)
hm3 = np.mean(attributions4[0][kmeans1.labels_==1,:,:]-attributions4B[0][kmeans1.labels_==1,:,:],axis=0)
hm4 = np.mean(attributions4[0][kmeans1.labels_==2,:,:]-attributions4B[0][kmeans1.labels_==2,:,:],axis=0)
hm5 = np.mean(attributions4[0][kmeans1.labels_==3,:,:]-attributions4B[0][kmeans1.labels_==3,:,:],axis=0)
hm6 = np.mean(attributions4[0][kmeans1.labels_==4,:,:]-attributions4B[0][kmeans1.labels_==4,:,:],axis=0)

kmeans2 = KMeans(n_clusters=5, random_state=0).fit(
                                                   np.reshape(attributions5[0],(np.shape(attributions5[0])[0],np.shape(attributions5[0])[1]*np.shape(attributions5[0])[2]))-np.reshape(attributions5B[0],(np.shape(attributions5B[0])[0],np.shape(attributions5B[0])[1]*np.shape(attributions5B[0])[2])) )
kmeans2.labels_

hm7 = np.mean(attributions5[0]-attributions5B[0],axis=0)
hm8 = np.mean(attributions5[0][kmeans2.labels_==0,:,:]-attributions5B[0][kmeans2.labels_==0,:,:],axis=0)
hm9 = np.mean(attributions5[0][kmeans2.labels_==1,:,:]-attributions5B[0][kmeans2.labels_==1,:,:],axis=0)
hm10 = np.mean(attributions5[0][kmeans2.labels_==2,:,:]-attributions5B[0][kmeans2.labels_==2,:,:],axis=0)
hm11 = np.mean(attributions5[0][kmeans2.labels_==3,:,:]-attributions5B[0][kmeans2.labels_==3,:,:],axis=0)
hm12 = np.mean(attributions5[0][kmeans2.labels_==4,:,:]-attributions5B[0][kmeans2.labels_==4,:,:],axis=0)


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



kmeans3 = KMeans(n_clusters=5, random_state=0).fit(
                                                   np.reshape(attributions4[0],(np.shape(attributions4[0])[0],np.shape(attributions4[0])[1]*np.shape(attributions4[0])[2]))- np.reshape(attributions4B[0],(np.shape(attributions4B[0])[0],np.shape(attributions4B[0])[1]*np.shape(attributions4B[0])[2])) )
kmeans3.labels_

hm13 = np.mean(attributions4[0]-attributions4B[0],axis=0)
hm14 = np.mean(attributions4[0][kmeans3.labels_==0,:,:]-attributions4B[0][kmeans3.labels_==0,:,:],axis=0)
hm15 = np.mean(attributions4[0][kmeans3.labels_==1,:,:]-attributions4B[0][kmeans3.labels_==1,:,:],axis=0)
hm16 = np.mean(attributions4[0][kmeans3.labels_==2,:,:]-attributions4B[0][kmeans3.labels_==2,:,:],axis=0)
hm17 = np.mean(attributions4[0][kmeans3.labels_==3,:,:]-attributions4B[0][kmeans3.labels_==3,:,:],axis=0)
hm18 = np.mean(attributions4[0][kmeans3.labels_==4,:,:]-attributions4B[0][kmeans3.labels_==4,:,:],axis=0)


kmeans4 = KMeans(n_clusters=5, random_state=0).fit(
                                                   np.reshape(attributions5[0],(np.shape(attributions5[0])[0],np.shape(attributions5[0])[1]*np.shape(attributions5[0])[2]))-np.reshape(attributions5B[0],(np.shape(attributions5B[0])[0],np.shape(attributions5B[0])[1]*np.shape(attributions5B[0])[2])) )
kmeans4.labels_

hm19 = np.mean(attributions5[0]-attributions5B[0],axis=0)
hm20 = np.mean(attributions5[0][kmeans4.labels_==0,:,:]-attributions5B[0][kmeans4.labels_==0,:,:],axis=0)
hm21 = np.mean(attributions5[0][kmeans4.labels_==1,:,:]-attributions5B[0][kmeans4.labels_==1,:,:],axis=0)
hm22 = np.mean(attributions5[0][kmeans4.labels_==2,:,:]-attributions5B[0][kmeans4.labels_==2,:,:],axis=0)
hm23 = np.mean(attributions5[0][kmeans4.labels_==3,:,:]-attributions5B[0][kmeans4.labels_==3,:,:],axis=0)
hm24 = np.mean(attributions5[0][kmeans4.labels_==4,:,:]-attributions5B[0][kmeans4.labels_==4,:,:],axis=0)

labs=["Ecto_H2A.ub","Ecto_H2A.x.3","Ecto_H3K27ac","Ecto_H3K27me3","Ecto_H3K36me3","Ecto_H3K4me1","Ecto_H3K4me2","Ecto_H3K4me3","Ecto_H3K79me3","Ecto_H3K9me3","Endo_H2A.ub","Endo_H2A.x.3","Endo_H3K27ac","Endo_H3K27me3","Endo_H3K36me3","Endo_H3K4me1","Endo_H3K4me2","Endo_H3K4me3","Endo_H3K79me3","Endo_H3K9me3","GC.bw.tab","Methylation"]


plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm1), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.0,0.001)
plt.colorbar()
plt.savefig('EndoChr/OnTP_Deltac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm2), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0,0.001)
plt.colorbar()
plt.savefig('EndoChr/ONTP_Deltac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm3), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0,0.001)
plt.colorbar()
plt.savefig('EndoChr/ONTP_Deltac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm4), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0,0.001)
plt.colorbar()
plt.savefig('EndoChr/ONTP_Deltac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm5), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0,0.001)
plt.colorbar()
plt.savefig('EndoChr/ONTP_Deltac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm6), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.0,0.001)
plt.colorbar()
plt.savefig('EndoChr/ONTP_Deltac_cl5_deeper_seed=2.pdf')
plt.close()


plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm7), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_Deltac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm8), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_Deltac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm9), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_Deltac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm10), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_Deltac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm11), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_Deltac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm12), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/OffTP_Deltac_cl5_deeper_seed=2.pdf')
plt.close()










plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm13), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.002)
plt.colorbar()
plt.savefig('EndoChr/RUTP_Deltac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm14), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.002)
plt.colorbar()
plt.savefig('EndoChr/RUTP_Deltac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm15), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.002)
plt.colorbar()
plt.savefig('EndoChr/RUTP_Deltac_cl2_deeper_seed=2.pdf')
plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm16), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.002)
plt.colorbar()
plt.savefig('EndoChr/RUTP_Deltac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm17), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.002)
plt.colorbar()
plt.savefig('EndoChr/RUTP_Deltac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm18), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.002)
plt.colorbar()
plt.savefig('EndoChr/RUTP_Deltac_cl5_deeper_seed=2.pdf')
plt.close()



plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm19), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/RDTP_Deltac_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm20), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/RDTP_Deltac_cl1_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm21), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/RDTP_Deltac_cl2_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm22), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/RDTP_Deltac_cl3_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm23), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/RDTP_Deltac_cl4_deeper_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm24), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChr/RDTP_Deltac_cl5_deeper_seed=2.pdf')
plt.close()

np.save('EndoChr/hm1_delta_seed=2.npy', hm1)
np.save('EndoChr/hm2_delta_seed=2.npy', hm2)
np.save('EndoChr/hm3_delta_seed=2.npy', hm3)
np.save('EndoChr/hm4_delta_seed=2.npy', hm4)
np.save('EndoChr/hm5_delta_seed=2.npy', hm4)
np.save('EndoChr/hm6_delta_seed=2.npy', hm6)
np.save('EndoChr/hm7_delta_seed=2.npy', hm7)
np.save('EndoChr/hm8_delta_seed=2.npy', hm8)
np.save('EndoChr/hm9_delta_seed=2.npy', hm9)
np.save('EndoChr/hm10_delta_seed=2.npy', hm10)
np.save('EndoChr/hm11_delta_seed=2.npy', hm11)
np.save('EndoChr/hm12_delta_seed=2.npy', hm12)
np.save('EndoChr/hm13_delta_seed=2.npy', hm13)
np.save('EndoChr/hm14_delta_seed=2.npy', hm14)
np.save('EndoChr/hm15_delta_seed=2.npy', hm15)
np.save('EndoChr/hm16_delta_seed=2.npy', hm16)
np.save('EndoChr/hm17_delta_seed=2.npy', hm17)
np.save('EndoChr/hm18_delta_seed=2.npy', hm18)
np.save('EndoChr/hm19_delta_seed=2.npy', hm19)
np.save('EndoChr/hm20_delta_seed=2.npy', hm20)
np.save('EndoChr/hm21_delta_seed=2.npy', hm21)
np.save('EndoChr/hm22_delta_seed=2.npy', hm22)
np.save('EndoChr/hm23_delta_seed=2.npy', hm23)
np.save('EndoChr/hm24_delta_seed=2.npy', hm24)


dsadsadasdsad

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



