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


Output2 = genfromtxt('BW_SomiteEndo_All/allpeaks_labelled_cum_final.se.bed',delimiter='\t',dtype=None)
#EndoNT_On_all.inde.bed',delimiter='\t',dtype=None)

#Load in the list of genes
onlist = genfromtxt('ChIP_SomiteEndo_All/OnMemFDRtGenes.txt',delimiter='\t',dtype=None)
downlist = genfromtxt('ChIP_SomiteEndo_All/ReprogrammedDowntGenes.txt',delimiter='\t',dtype=None)
offlist = genfromtxt('ChIP_SomiteEndo_All/OffMemFDRtGenes.txt',delimiter='\t',dtype=None)
uplist = genfromtxt('ChIP_SomiteEndo_All/ReprogrammedUptGenes.txt',delimiter='\t',dtype=None)
complist = genfromtxt('ChIP_SomiteEndo_All/ComplementarySettGenes.txt',delimiter='\t',dtype=None)



#Generate the Y matrix
Yalt = np.zeros((np.shape(Output2)[0],5))
Xexpa = np.zeros((np.shape(Output2)[0],1,2))
Xrawexp = np.zeros((np.shape(Output2)[0],7))

for i in range(0, np.shape(Yalt)[0]):
      Yalt[i,0] = (np.intersect1d(Output2['f3'][i],onlist)).size
      Yalt[i,1] = (np.intersect1d(Output2['f3'][i],downlist)).size
      Yalt[i,2] = (np.intersect1d(Output2['f3'][i],offlist)).size
      Yalt[i,3] = (np.intersect1d(Output2['f3'][i],uplist)).size
      Yalt[i,4] = 1-max(Yalt[i,0:4]) 

explist = genfromtxt('ChIP_SomiteEndo_All/edger_de_pairing_BIGtable_endoIVF_somiteDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)

explist1 = genfromtxt('ChIP_SomiteEndo_All/edger_de_pairing_BIGtable_endoIVF_somiteDonor_endoNT_plus_DE.nan.csv',delimiter='\t',dtype=None)

#explist = genfromtxt('ChIP_SomiteEndo_All/edger_de_pairing_BIGtable_ectoIVF_somiteDonor_ectoNT_plus_DE.p.csv',delimiter='\t',dtype=None)
#explist = genfromtxt('ChIP_SomiteEcto_All/edger_de_pairing_BIGtable_endoIVF_ectoDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)

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

fil=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All/*bwe.tab'))

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

model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAll/Model_deeper_seed=2.h5')
model1 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAll/Model_deeper_neg_seed=2.h5')
model2 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAll/Model_deeper_pos_seed=2.h5')


#model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
#model1 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_neg.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
#model2 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_pos.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)

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

np.save('EndoFullSomEndo/OnTP_seed=2.npy', OnTP)
np.save('EndoFullSomEndo/OnFP_seed=2.npy', OnFP)
np.save('EndoFullSomEndo/OnTN_seed=2.npy', OnTN)
np.save('EndoFullSomEndo/OnFN_seed=2.npy', OnFN)
np.save('EndoFullSomEndo/OffTP_seed=2.npy', OffTP)
np.save('EndoFullSomEndo/OffFP_seed=2.npy', OffFP)
np.save('EndoFullSomEndo/OffTN_seed=2.npy', OffTN)
np.save('EndoFullSomEndo/OffFN_seed=2.npy', OffFN)
np.save('EndoFullSomEndo/RUTP_seed=2.npy', RUTP)
np.save('EndoFullSomEndo/RUFP_seed=2.npy', RUFP)
np.save('EndoFullSomEndo/RUTN_seed=2.npy', RUTN)
np.save('EndoFullSomEndo/RUFN_seed=2.npy', RUFN)
np.save('EndoFullSomEndo/RDTP_seed=2.npy', RDTP)
np.save('EndoFullSomEndo/RDFP_seed=2.npy', RDFP)
np.save('EndoFullSomEndo/RDTN_seed=2.npy', RDTN)
np.save('EndoFullSomEndo/RDFN_seed=2.npy', RDFN)


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

np.save('EndoFullSomEndo/On_attributions_OnTP_seed=2.npy', attributions4)
np.save('EndoFullSomEndo/Down_attributions_OnTP_seed=2.npy', attributions4B)
np.save('EndoFullSomEndo/Off_attributions_OffTP_seed=2.npy', attributions5)
np.save('EndoFullSomEndo/Up_attributions_OffTP_seed=2.npy', attributions5B)


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

np.save('EndoFullSomEndo/Down_attributions_RUTP_seed=2.npy', attributions4)
np.save('EndoFullSomEndo/On_attributions_RUTP_seed=2.npy', attributions4B)
np.save('EndoFullSomEndo/Up_attributions_RDTP_seed=2.npy', attributions5)
np.save('EndoFullSomEndo/Off_attributions_RDTP_seed=2.npy', attributions5B)

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
#import matplotlib.pyplot as plt
#from matplotlib import colors as mcolors
#plt.switch_backend('agg')

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm1), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OnTP_ac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm2), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/ONTP_ac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm3), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/ONTP_ac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm4), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/ONTP_ac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm5), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/ONTP_ac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm6), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/ONTP_ac_cl5_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm7), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.006,0.006)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_ac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm8), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.006,0.006)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_ac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm9), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.006,0.006)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_ac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm10), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.006,0.006)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_ac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm11), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.006,0.006)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_ac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm12), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.006,0.006)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_ac_cl5_seed=2.pdf')
plt.close()


plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm13), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_ac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm14), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_ac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm15), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_ac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm16), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_ac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm17), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_ac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm18), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_ac_cl5_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm19), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_ac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm20), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_ac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm21), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_ac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm22), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_ac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm23), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_ac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm24), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.003)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_ac_cl5_seed=2.pdf')
plt.close()



np.save('EndoFullSomEndo/hm1_seed=2.npy', hm1)
np.save('EndoFullSomEndo/hm2_seed=2.npy', hm2)
np.save('EndoFullSomEndo/hm3_seed=2.npy', hm3)
np.save('EndoFullSomEndo/hm4_seed=2.npy', hm4)
np.save('EndoFullSomEndo/hm5_seed=2.npy', hm5)
np.save('EndoFullSomEndo/hm6_seed=2.npy', hm6)
np.save('EndoFullSomEndo/hm7_seed=2.npy', hm7)
np.save('EndoFullSomEndo/hm8_seed=2.npy', hm8)
np.save('EndoFullSomEndo/hm9_seed=2.npy', hm9)
np.save('EndoFullSomEndo/hm10_seed=2.npy', hm10)
np.save('EndoFullSomEndo/hm11_seed=2.npy', hm11)
np.save('EndoFullSomEndo/hm12_seed=2.npy', hm12)
np.save('EndoFullSomEndo/hm13_seed=2.npy', hm13)
np.save('EndoFullSomEndo/hm14_seed=2.npy', hm14)
np.save('EndoFullSomEndo/hm15_seed=2.npy', hm15)
np.save('EndoFullSomEndo/hm16_seed=2.npy', hm16)
np.save('EndoFullSomEndo/hm17_seed=2.npy', hm17)
np.save('EndoFullSomEndo/hm18_seed=2.npy', hm18)
np.save('EndoFullSomEndo/hm19_seed=2.npy', hm19)
np.save('EndoFullSomEndo/hm20_seed=2.npy', hm20)
np.save('EndoFullSomEndo/hm21_seed=2.npy', hm21)
np.save('EndoFullSomEndo/hm22_seed=2.npy', hm22)
np.save('EndoFullSomEndo/hm23_seed=2.npy', hm23)
np.save('EndoFullSomEndo/hm24_seed=2.npy', hm24)



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
plt.clim(-0.003,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OnTP_Deltac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm2), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
#plt.clim(-0.001,0.003)
plt.clim(-0.003,0.004)

plt.colorbar()
plt.savefig('EndoFullSomEndo/ONTP_Deltac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm3), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
#plt.clim(-0.001,0.003)
plt.clim(-0.003,0.004)

plt.colorbar()
plt.savefig('EndoFullSomEndo/ONTP_Deltac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm4), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
#plt.clim(-0.001,0.003)
plt.clim(-0.003,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/ONTP_Deltac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm5), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.003,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/ONTP_Deltac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm6), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.003,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/ONTP_Deltac_cl5_seed=2.pdf')
plt.close()


plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm7), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.006,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_Deltac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm8), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.006,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_Deltac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm9), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.006,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_Deltac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm10), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.006,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_Deltac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm11), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.006,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_Deltac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm12), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.006,0.004)
plt.colorbar()
plt.savefig('EndoFullSomEndo/OffTP_Deltac_cl5_seed=2.pdf')
plt.close()










plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm13), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.004,0.001)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_Deltac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm14), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.001)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_Deltac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm15), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.001)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_Deltac_cl2_seed=2.pdf')
plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm16), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.001)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_Deltac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm17), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.001)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_Deltac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm18), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.001)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RUTP_Deltac_cl5_seed=2.pdf')
plt.close()



plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm19), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_Deltac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm20), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_Deltac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm21), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_Deltac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm22), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)

plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_Deltac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm23), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_Deltac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm24), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.002)
plt.colorbar()
plt.savefig('EndoFullSomEndo/RDTP_Deltac_cl5_seed=2.pdf')
plt.close()



np.save('EndoFullSomEndo/hm1_delta_seed=2.npy', hm1)
np.save('EndoFullSomEndo/hm2_delta_seed=2.npy', hm2)
np.save('EndoFullSomEndo/hm3_delta_seed=2.npy', hm3)
np.save('EndoFullSomEndo/hm4_delta_seed=2.npy', hm4)
np.save('EndoFullSomEndo/hm5_delta_seed=2.npy', hm5)
np.save('EndoFullSomEndo/hm6_delta_seed=2.npy', hm6)
np.save('EndoFullSomEndo/hm7_delta_seed=2.npy', hm7)
np.save('EndoFullSomEndo/hm8_delta_seed=2.npy', hm8)
np.save('EndoFullSomEndo/hm9_delta_seed=2.npy', hm9)
np.save('EndoFullSomEndo/hm10_delta_seed=2.npy', hm10)
np.save('EndoFullSomEndo/hm11_delta_seed=2.npy', hm11)
np.save('EndoFullSomEndo/hm12_delta_seed=2.npy', hm12)
np.save('EndoFullSomEndo/hm13_delta_seed=2.npy', hm13)
np.save('EndoFullSomEndo/hm14_delta_seed=2.npy', hm14)
np.save('EndoFullSomEndo/hm15_delta_seed=2.npy', hm15)
np.save('EndoFullSomEndo/hm16_delta_seed=2.npy', hm16)
np.save('EndoFullSomEndo/hm17_delta_seed=2.npy', hm17)
np.save('EndoFullSomEndo/hm18_delta_seed=2.npy', hm18)
np.save('EndoFullSomEndo/hm19_delta_seed=2.npy', hm19)
np.save('EndoFullSomEndo/hm20_delta_seed=2.npy', hm20)
np.save('EndoFullSomEndo/hm21_delta_seed=2.npy', hm21)
np.save('EndoFullSomEndo/hm22_delta_seed=2.npy', hm22)
np.save('EndoFullSomEndo/hm23_delta_seed=2.npy', hm23)
np.save('EndoFullSomEndo/hm24_delta_seed=2.npy', hm24)

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


np.save('EndoChrSomEndo/OnTP_seed=2.npy', OnTP)
np.save('EndoChrSomEndo/OnFP_seed=2.npy', OnFP)
np.save('EndoChrSomEndo/OnTN_seed=2.npy', OnTN)
np.save('EndoChrSomEndo/OnFN_seed=2.npy', OnFN)
np.save('EndoChrSomEndo/OffTP_seed=2.npy', OffTP)
np.save('EndoChrSomEndo/OffFP_seed=2.npy', OffFP)
np.save('EndoChrSomEndo/OffTN_seed=2.npy', OffTN)
np.save('EndoChrSomEndo/OffFN_seed=2.npy', OffFN)
np.save('EndoChrSomEndo/RUTP_seed=2.npy', RUTP)
np.save('EndoChrSomEndo/RUFP_seed=2.npy', RUFP)
np.save('EndoChrSomEndo/RUTN_seed=2.npy', RUTN)
np.save('EndoChrSomEndo/RUFN_seed=2.npy', RUFN)
np.save('EndoChrSomEndo/RDTP_seed=2.npy', RDTP)
np.save('EndoChrSomEndo/RDFP_seed=2.npy', RDFP)
np.save('EndoChrSomEndo/RDTN_seed=2.npy', RDTN)
np.save('EndoChrSomEndo/RDFN_seed=2.npy', RDFN)



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

np.save('EndoChrSomEndo/On_attributions_OnTP_seed=2.npy', attributions4)
np.save('EndoChrSomEndo/Down_attributions_OnTP_seed=2.npy', attributions4B)
np.save('EndoChrSomEndo/Off_attributions_OffTP_seed=2.npy', attributions5)
np.save('EndoChrSomEndo/Up_attributions_OffTP_seed=2.npy', attributions5B)



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


np.save('EndoChrSomEndo/Down_attributions_RUTP_seed=2.npy', attributions4)
np.save('EndoChrSomEndo/On_attributions_RUTP_seed=2.npy', attributions4B)
np.save('EndoChrSomEndo/Up_attributions_RDTP_seed=2.npy', attributions5)
np.save('EndoChrSomEndo/Off_attributions_RDTP_seed=2.npy', attributions5B)


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
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OnTP_ac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm2), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.003)
#plt.clim(-0.003,0.001)
plt.colorbar()
plt.savefig('EndoChrSomEndo/ONTP_ac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm3), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/ONTP_ac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm4), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/ONTP_ac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm5), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/ONTP_ac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm6), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/ONTP_ac_cl5_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm7), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.004,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_ac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm8), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_ac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm9), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_ac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm10), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_ac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm11), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_ac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm12), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_ac_cl5_seed=2.pdf')
plt.close()




plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm13), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_ac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm14), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_ac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm15), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_ac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm16), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_ac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm17), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_ac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm18), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_ac_cl5_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm19), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RDTP_ac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm20), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RDTP_ac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm21), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RDTP_ac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm22), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RDTP_ac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm23), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RDTP_ac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm24), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RDTP_ac_cl5_seed=2.pdf')
plt.close()



np.save('EndoChrSomEndo/hm1_seed=2.npy', hm1)
np.save('EndoChrSomEndo/hm2_seed=2.npy', hm2)
np.save('EndoChrSomEndo/hm3_seed=2.npy', hm3)
np.save('EndoChrSomEndo/hm4_seed=2.npy', hm4)
np.save('EndoChrSomEndo/hm5_seed=2.npy', hm4)
np.save('EndoChrSomEndo/hm6_seed=2.npy', hm6)
np.save('EndoChrSomEndo/hm7_seed=2.npy', hm7)
np.save('EndoChrSomEndo/hm8_seed=2.npy', hm8)
np.save('EndoChrSomEndo/hm9_seed=2.npy', hm9)
np.save('EndoChrSomEndo/hm10_seed=2.npy', hm10)
np.save('EndoChrSomEndo/hm11_seed=2.npy', hm11)
np.save('EndoChrSomEndo/hm12_seed=2.npy', hm12)
np.save('EndoChrSomEndo/hm13_seed=2.npy', hm13)
np.save('EndoChrSomEndo/hm14_seed=2.npy', hm14)
np.save('EndoChrSomEndo/hm15_seed=2.npy', hm15)
np.save('EndoChrSomEndo/hm16_seed=2.npy', hm16)
np.save('EndoChrSomEndo/hm17_seed=2.npy', hm17)
np.save('EndoChrSomEndo/hm18_seed=2.npy', hm18)
np.save('EndoChrSomEndo/hm19_seed=2.npy', hm19)
np.save('EndoChrSomEndo/hm20_seed=2.npy', hm20)
np.save('EndoChrSomEndo/hm21_seed=2.npy', hm21)
np.save('EndoChrSomEndo/hm22_seed=2.npy', hm22)
np.save('EndoChrSomEndo/hm23_seed=2.npy', hm23)
np.save('EndoChrSomEndo/hm24_seed=2.npy', hm24)



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
plt.clim(-0.003,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OnTP_Deltac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm2), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.003,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/ONTP_Deltac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm3), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.003,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/ONTP_Deltac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm4), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.003,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/ONTP_Deltac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm5), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.003,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/ONTP_Deltac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm6), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.003,0.003)
plt.colorbar()
plt.savefig('EndoChrSomEndo/ONTP_Deltac_cl5_seed=2.pdf')
plt.close()


plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm7), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_Deltac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm8), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_Deltac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm9), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_Deltac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm10), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_Deltac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm11), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_Deltac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm12), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.004,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/OffTP_Deltac_cl5_seed=2.pdf')
plt.close()










plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm13), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.001,0.001)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_Deltac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm14), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.001)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_Deltac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm15), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.001)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_Deltac_cl2_seed=2.pdf')
plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm16), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.001)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_Deltac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm17), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.001)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_Deltac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm18), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.001,0.001)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RUTP_Deltac_cl5_seed=2.pdf')
plt.close()



plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm19), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
plt.clim(-0.002,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RDTP_Deltac_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm20), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.004)
plt.colorbar()
plt.savefig('EndoHChrSomEndo/RDTP_Deltac_cl1_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm21), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RDTP_Deltac_cl2_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm22), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RDTP_Deltac_cl3_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm23), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.004)
plt.colorbar()
plt.savefig('EndoChrSomEndo/RDTP_Deltac_cl4_seed=2.pdf')
plt.close()

plt.figure()
plt.subplot(1, 1, 1)
plt.imshow(np.transpose(hm24), interpolation='nearest', aspect='auto')
plt.yticks(np.arange(22), labs, fontsize=6,  rotation=20)
#plt.clim(-0.003,0.001)
plt.clim(-0.002,0.004)

plt.colorbar()
plt.savefig('EndoChrSomEndo/RDTP_Deltac_cl5_seed=2.pdf')
plt.close()




np.save('EndoChrSomEndo/hm1_delta_seed=2.npy', hm1)
np.save('EndoChrSomEndo/hm2_delta_seed=2.npy', hm2)
np.save('EndoChrSomEndo/hm3_delta_seed=2.npy', hm3)
np.save('EndoChrSomEndo/hm4_delta_seed=2.npy', hm4)
np.save('EndoChrSomEndo/hm5_delta_seed=2.npy', hm4)
np.save('EndoChrSomEndo/hm6_delta_seed=2.npy', hm6)
np.save('EndoChrSomEndo/hm7_delta_seed=2.npy', hm7)
np.save('EndoChrSomEndo/hm8_delta_seed=2.npy', hm8)
np.save('EndoChrSomEndo/hm9_delta_seed=2.npy', hm9)
np.save('EndoChrSomEndo/hm10_delta_seed=2.npy', hm10)
np.save('EndoChrSomEndo/hm11_delta_seed=2.npy', hm11)
np.save('EndoChrSomEndo/hm12_delta_seed=2.npy', hm12)
np.save('EndoChrSomEndo/hm13_delta_seed=2.npy', hm13)
np.save('EndoChrSomEndo/hm14_delta_seed=2.npy', hm14)
np.save('EndoChrSomEndo/hm15_delta_seed=2.npy', hm15)
np.save('EndoChrSomEndo/hm16_delta_seed=2.npy', hm16)
np.save('EndoChrSomEndo/hm17_delta_seed=2.npy', hm17)
np.save('EndoChrSomEndo/hm18_delta_seed=2.npy', hm18)
np.save('EndoChrSomEndo/hm19_delta_seed=2.npy', hm19)
np.save('EndoChrSomEndo/hm20_delta_seed=2.npy', hm20)
np.save('EndoChrSomEndo/hm21_delta_seed=2.npy', hm21)
np.save('EndoChrSomEndo/hm22_delta_seed=2.npy', hm22)
np.save('EndoChrSomEndo/hm23_delta_seed=2.npy', hm23)
np.save('EndoChrSomEndo/hm24_delta_seed=2.npy', hm24)




















