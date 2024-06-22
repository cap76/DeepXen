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

#List of genes locations
Output2 = genfromtxt('BW_Endo_All2/allpeaks_labelled_cum_final.se.bed',delimiter='\t',dtype=None)

#Load in the list of genes
onlist = genfromtxt('ChIP_Endo_All/OnMemFDRtGenes.txt',delimiter='\t',dtype=None)
downlist = genfromtxt('ChIP_Endo_All/ReprogrammedDowntGenes.txt',delimiter='\t',dtype=None)
offlist = genfromtxt('ChIP_Endo_All/OffMemFDRtGenes.txt',delimiter='\t',dtype=None)
uplist = genfromtxt('ChIP_Endo_All/ReprogrammedUptGenes.txt',delimiter='\t',dtype=None)
complist = genfromtxt('ChIP_Endo_All/ComplementarySettGenes.txt',delimiter='\t',dtype=None)

#foxlist = genfromtxt('BW_Endo_All2/FOXA1ChIP.txt',delimiter='\t',dtype=None)
#enhlist = genfromtxt('BW_Endo_All2/enhancers.bed',delimiter='\t',dtype=None)


#Generate the Y matrix
Yalt = np.zeros((np.shape(Output2)[0],5))
Xexpa = np.zeros((np.shape(Output2)[0],1,2))

for i in range(0, np.shape(Yalt)[0]):
      Yalt[i,0] = (np.intersect1d(Output2['f3'][i],onlist)).size
      Yalt[i,1] = (np.intersect1d(Output2['f3'][i],downlist)).size
      Yalt[i,2] = (np.intersect1d(Output2['f3'][i],offlist)).size
      Yalt[i,3] = (np.intersect1d(Output2['f3'][i],uplist)).size
      Yalt[i,4] = 1-max(Yalt[i,0:4]) 

explist = genfromtxt('ChIP_Endo_All/edger_de_pairing_BIGtable_endoIVF_ectoDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)

for i in range(0, np.shape(Yalt)[0]):
   arr_ind = np.where(Output2['f3'][i] == explist['f0'])
   Xexpa[i,0,0] = np.log2(explist['f1'][arr_ind].astype('f')+1)
   Xexpa[i,0,1] = np.log2(explist['f2'][arr_ind].astype('f')+1)


np.nan_to_num(Xexpa,copy=False)

#List the files
fila=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/*Input*bwe.tab'))
filb=glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/*GC*bwe.tab')
filc=glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/*Methylation*bwe.tab')
fil = fila + filb + filc
#Load in all the expression data
encoder = LabelEncoder()
meh = encoder.fit(Output2['f0'])
Xchr = np.zeros((np.shape(Output2)[0],600,np.shape(fil)[0]))
for k in range(0, np.shape(fil)[0]):
        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
        Xchr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output

#Get the TFs
fil2=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/Motif*top*.tab'))
Outputlab = genfromtxt('BW_Endo_All2/Motif_ap1.bed.bot.bw_oute.bed',delimiter='\t',dtype=None)
XTF1 = np.zeros((np.shape(Output2)[0],600,np.shape(fil2)[0]))
XTF2 = np.zeros((np.shape(Output2)[0],600,np.shape(fil2)[0]))

#Chromatin mods
for k in range(0, np.shape(fil)[0]):
        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
        Xchr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output

#Better to frame this as an outer join? For now just get indexes of where the TF region maps to the peak region
indexss = np.zeros((np.shape(Outputlab)[0],),dtype=np.int8)
for k in range(0, np.shape(Outputlab)[0]):
       indexss[k] = np.where((Outputlab['f0'][k]==Output2['f0']) & (Outputlab['f1'][k]==Output2['f1']) & (Outputlab['f2'][k]==Output2['f2']))[0][0]

#Load in TF regions
for k in range(0, np.shape(fil2)[0]):
      try:
         Output = genfromtxt(fil2[k],delimiter='\t',dtype='f', skip_header=3)
         XTF1[indexss.astype(np.int64),:,0] = Output
         Output = genfromtxt(fil2[k].replace("top","bot"),delimiter='\t',dtype='f', skip_header=3)
         XTF2[indexss.astype(np.int64),:,0] = Output
      except:
         print(k)
         print(fil2[k])
         continue

np.nan_to_num(XTF2,copy=False)
np.nan_to_num(XTF1,copy=False)

#Now merge the three datasets?

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

visible1 = Input(shape=(600,np.shape(fil)[0]))
#Convolution layers
conv3_1 = Conv1D(64, kernel_size=16, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(visible1)
pool3_1 = MaxPooling1D(pool_size=4)(conv3_1)
conv3_2 = Conv1D(64, kernel_size=16, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool3_1)
pool3_2 = MaxPooling1D(pool_size=4)(conv3_2)
conv3_3 = Conv1D(64, kernel_size=16, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool3_2)
pool3_3 = MaxPooling1D(pool_size=4)(conv3_3)

layer3_1b = Dropout(0.4)(pool3_3)
flat3 = Flatten()(layer3_1b)

conv4_1 = Conv1D(64, kernel_size=32, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(visible1)
pool4_1 = MaxPooling1D(pool_size=4)(conv4_1)
conv4_2 = Conv1D(64, kernel_size=32, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool4_1)
pool4_2 = MaxPooling1D(pool_size=4)(conv4_2)
conv4_3 = Conv1D(64, kernel_size=32, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool4_2)
pool4_3 = MaxPooling1D(pool_size=4)(conv4_3)

layer4_1b = Dropout(0.4)(pool4_3)
flat4 = Flatten()(layer4_1b)


#Global scan
conv6_1 = Conv1D(64, kernel_size=400, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(visible1)
pool6_1 = MaxPooling1D(pool_size=4)(conv6_1)
layer6_1b = Dropout(0.4)(pool6_1)
flat6 = Flatten()(layer6_1b)






visible1a = Input(shape=(600,np.shape(fil2)[0]))
#Convolution layers
conv3_1a = Conv1D(64, kernel_size=16, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(visible1a)
pool3_1a = MaxPooling1D(pool_size=4)(conv3_1a)
conv3_2a = Conv1D(64, kernel_size=16, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool3_1a)
pool3_2a = MaxPooling1D(pool_size=4)(conv3_2a)
conv3_3a = Conv1D(64, kernel_size=16, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool3_2a)
pool3_3a = MaxPooling1D(pool_size=4)(conv3_3a)

layer3_1ba = Dropout(0.4)(pool3_3a)
flat3a = Flatten()(layer3_1ba)

conv4_1a = Conv1D(64, kernel_size=32, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(visible1a)
pool4_1a = MaxPooling1D(pool_size=4)(conv4_1a)
conv4_2a = Conv1D(64, kernel_size=32, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool4_1a)
pool4_2a = MaxPooling1D(pool_size=4)(conv4_2a)
conv4_3a = Conv1D(64, kernel_size=32, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool4_2a)
pool4_3a = MaxPooling1D(pool_size=4)(conv4_3a)

layer4_1ba = Dropout(0.4)(pool4_3a)
flat4a = Flatten()(layer4_1ba)


#Global scan
conv6_1a = Conv1D(64, kernel_size=400, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(visible1a)
pool6_1a = MaxPooling1D(pool_size=4)(conv6_1a)
layer6_1ba = Dropout(0.4)(pool6_1a)
flat6a = Flatten()(layer6_1ba)




visible1b = Input(shape=(600,np.shape(fil2)[0]))
#Convolution layers
conv3_1b = Conv1D(64, kernel_size=16, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(visible1b)
pool3_1b = MaxPooling1D(pool_size=4)(conv3_1b)
conv3_2b = Conv1D(64, kernel_size=16, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool3_1b)
pool3_2b = MaxPooling1D(pool_size=4)(conv3_2b)
conv3_3b = Conv1D(64, kernel_size=16, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool3_2b)
pool3_3b = MaxPooling1D(pool_size=4)(conv3_3b)

layer3_1bb = Dropout(0.4)(pool3_3b)
flat3b = Flatten()(layer3_1bb)

conv4_1b = Conv1D(64, kernel_size=32, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(visible1b)
pool4_1b = MaxPooling1D(pool_size=4)(conv4_1b)
conv4_2b = Conv1D(64, kernel_size=32, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool4_1b)
pool4_2b = MaxPooling1D(pool_size=4)(conv4_2b)
conv4_3b = Conv1D(64, kernel_size=32, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(pool4_2b)
pool4_3b = MaxPooling1D(pool_size=4)(conv4_3b)

layer4_1bb = Dropout(0.4)(pool4_3b)
flat4b = Flatten()(layer4_1bb)


#Global scan
conv6_1b = Conv1D(64, kernel_size=400, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True,padding='same')(visible1b)
pool6_1b = MaxPooling1D(pool_size=4)(conv6_1b)
layer6_1bb = Dropout(0.4)(pool6_1b)
flat6b = Flatten()(layer6_1bb)




#Auxillary input (fully connected) for IVF and Donor expression (and other motifs)
visible2 = Input(shape=(1,2))

layer5_1 = Dense(10, activation='relu')(visible2)
flat5 = Flatten()(layer5_1)

merge = concatenate([flat3,flat4,flat6, flat3a,flat4a,flat6a, flat3b,flat4b,flat6b, flat5])

# interpretation model
hidden1 = Dense(50, activation='relu')(merge)
hidden2 = Dense(20, activation='relu')(hidden1)

output1 = Dense(np.shape(Yalt)[1])(hidden2)
output  = Activation('softmax')(output1)

model = Model(inputs=[visible1,visible1a,visible1b,visible2], outputs=output)
# summarize layers
print(model.summary())
# plot graph


nclasses = np.sum(Yalt[trainset,:],0)
cwe = np.ndarray.max(nclasses) / nclasses 
class_weight={0: cwe[0], 1: cwe[1], 2: cwe[2], 3: cwe[3], 4: cwe[4]} #, 5: cwe[5]} #, 6: cwe[6]}

sgd = optimizers.SGD(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)

model.compile(loss="categorical_crossentropy", optimizer = sgd, metrics=['categorical_accuracy'])
callbacks = [GetBest(monitor='val_categorical_accuracy', verbose=1, mode='max')]

model.fit([Xchr[trainset,:,:],XTF1[trainset,:,:],XTF2[trainset,:,:],Xexp[trainset,:,:]], Yalt[trainset,:], validation_data=([Xchr[testset,:,:],XTF1[testset,:,:],XTF2[testset,:,:],Xexp[testset,:,:]], Yalt[testset,:]), epochs=1000, batch_size=1000,class_weight = class_weight,callbacks=callbacks)

model.save('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_TF.h5')

predictions = model.predict([Xchr,XTF1,XTF2,Xexp], batch_size=1000)

Scores = np.zeros((3,1))
Scores[0,0] = model.evaluate([Xchr[trainset,:,:],XTF1[trainset,:,:],XTF2[trainset,:,:],Xexp[trainset,:,:]], Yalt[trainset,:], batch_size=32)[1]
Scores[1,0] = model.evaluate([Xchr[testset,:,:],XTF1[testset,:,:],XTF2[testset,:,:],Xexp[testset,:,:]], Yalt[testset,:], batch_size=32)[1]
Scores[2,0] = model.evaluate([Xchr[valset,:,:],XTF1[valset,:,:],XTF2[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]

pd.DataFrame(Scores, columns=['Accuracy']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/prediction_scores_deeper_TF.csv')

df1 = pd.DataFrame(Output2['f0'],columns=['chr'])
df2 = pd.DataFrame(Output2['f1'],columns=['start'])
df3 = pd.DataFrame(Output2['f2'],columns=['end'])
df4 = pd.DataFrame(Output2['f3'],columns=['C1'])
df5 = pd.DataFrame(Output2['f4'],columns=['C2'])
df6 = pd.DataFrame(Output2['f5'],columns=['C3'])
df7 = pd.DataFrame(Output2['f6'],columns=['C3'])
df8 = pd.DataFrame(Output2['f7'],columns=['C4'])
df9 = pd.DataFrame(Output2['f8'],columns=['C5'])
df10 = pd.DataFrame(Output2['f9'],columns=['C6'])
df11 = pd.DataFrame(Output2['f10'],columns=['C7'])
df12 = pd.DataFrame(Yalt[:,0],columns=['C8'])
df13 = pd.DataFrame(Yalt[:,1],columns=['C9'])
df14 = pd.DataFrame(Yalt[:,2],columns=['C10'])
df15 = pd.DataFrame(Yalt[:,3],columns=['C11'])
df16 = pd.DataFrame(Yalt[:,4],columns=['C12'])

df14b = pd.DataFrame(predictions,columns=['pC1','pC2','pC3','pC4','pC5'])

prediction = np.concatenate((df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14,df15,df16, df14b),1)
pd.DataFrame(prediction, columns=['Chr','start','end','Gene','IVF','Donor','NT','FC','pVal','FC','pVal','ON','RepDown','Off','RepUp','Neg','pC1','pC2','pC3','pC4','pC5']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/prediction_deeper_TF.csv')


