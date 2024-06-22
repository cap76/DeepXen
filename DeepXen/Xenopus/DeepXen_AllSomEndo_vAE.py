# Multiple Inputs
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import metrics


from tensorflow.keras.utils import plot_model

from tensorflow.keras.utils import plot_model
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import Flatten
from tensorflow.keras.layers import Activation
from tensorflow.keras.layers import Dropout
from tensorflow.keras import optimizers
from tensorflow.keras import regularizers
from tensorflow.keras.layers import Conv2D
from tensorflow.keras.layers import Conv1D
from tensorflow.keras.layers import MaxPooling2D
from tensorflow.keras.layers import MaxPooling1D
from tensorflow.keras.layers import concatenate

from tensorflow.keras.utils import to_categorical
#from tensorflow.keras.utils import np_utils
from sklearn.preprocessing import LabelEncoder
from numpy import genfromtxt
import numpy as np
import pandas as pd
import glob
import os
import sys
import copy

np.random.seed(0)

#Load in the raw data (labels are encoded in column 4) and one-hot encode.
Output2 = genfromtxt('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All/allpeaks_labelled_cum_final.s.bed',delimiter='\t',dtype=None)

#Load in the list of genes
onlist = genfromtxt('ChIP_SomiteEndo_All/OnMemFDRtGenes.txt',delimiter='\t',dtype=None)
downlist = genfromtxt('ChIP_SomiteEndo_All/ReprogrammedDowntGenes.txt',delimiter='\t',dtype=None)
offlist = genfromtxt('ChIP_SomiteEndo_All/OffMemFDRtGenes.txt',delimiter='\t',dtype=None)
uplist = genfromtxt('ChIP_SomiteEndo_All/ReprogrammedUptGenes.txt',delimiter='\t',dtype=None)
complist = genfromtxt('ChIP_SomiteEndo_All/ComplementarySettGenes.txt',delimiter='\t',dtype=None)

foxlist = genfromtxt('BW_Endo_All2/FOXA1ChIP.txt',delimiter='\t',dtype=None)

enhlist = genfromtxt('BW_Endo_All2/enhancers.bed',delimiter='\t',dtype=None)


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


explist = genfromtxt('ChIP_SomiteEndo_All/edger_de_pairing_BIGtable_endoIVF_somiteDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)

for i in range(0, np.shape(Yalt)[0]):
   arr_ind = np.where(Output2['f3'][i] == explist['f0'])
   Xexpa[i,0,0] = np.log2(explist['f1'][arr_ind].astype('f')+1)
   Xexpa[i,0,1] = np.log2(explist['f2'][arr_ind].astype('f')+1)

np.nan_to_num(Xexpa,copy=False)

fil=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All/*.tab'))

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
      Xchr[:,0:100,k] = ( Xchr[:,0:600,k] - np.nanmean(Xchr[trainset,0:600,k]) ) / np.nanstd(Xchr[trainset,0:600,k])

#Now do the thing

visible1 = Input(shape=(600,np.shape(fil)[0]))
#Convolution layers
conv3_1 = Conv1D(64, kernel_size=16, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True)(visible1)
pool3_1 = MaxPooling1D(pool_size=4)(conv3_1)
layer3_1b = Dropout(0.4)(pool3_1)
flat3 = Flatten()(layer3_1b)
conv4_1 = Conv1D(64, kernel_size=32, activation='relu',kernel_regularizer=regularizers.l1_l2(0.1),use_bias=True)(visible1)
pool4_1 = MaxPooling1D(pool_size=4)(conv4_1)
layer4_1b = Dropout(0.4)(pool4_1)
flat4 = Flatten()(layer4_1b)

#Auxillary input (fully connected) for IVF and Donor expression (and other motifs)
visible2 = Input(shape=(1,2))
layer5_1 = Dense(10, activation='relu')(visible2)
flat5 = Flatten()(layer5_1)
merge = concatenate([flat3,flat4,flat5])

# interpretation model
hidden1 = Dense(50, activation='relu')(merge)
hidden2 = Dense(20, activation='relu')(hidden1)
output1 = Dense(np.shape(Yalt)[1])(hidden2)
output  = Activation('softmax')(output1)

model = Model(inputs=[visible1,visible2], outputs=output)
# summarize layers
print(model.summary())
# plot graph


nclasses = np.sum(Yalt[trainset,:],0)
cwe = np.ndarray.max(nclasses) / nclasses 
class_weight={0: cwe[0], 1: cwe[1], 2: cwe[2], 3: cwe[3], 4: cwe[4]} #, 5: cwe[5]} #, 6: cwe[6]}

#sgd = optimizers.SGD(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)
#model.compile(loss="categorical_crossentropy", optimizer = sgd, metrics=['categorical_accuracy'])
#callbacks = [GetBest(monitor='val_categorical_accuracy', verbose=1, mode='max')]
#model.fit([Xchr[trainset,:,:],Xexp[trainset,:,:]], Yalt[trainset,:], validation_data=([Xchr[testset,:,:],Xexp[testset,:,:]], Yalt[testset,:]), epochs=1000, batch_size=1000,class_weight = class_weight,callbacks=callbacks)
#model.save('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomEndoAll/Model.h5')


#predictions = model.predict([Xchr,Xexp], batch_size=1000)
#Scores = np.zeros((3,1))
#Scores[0,0] = model.evaluate([Xchr[trainset,:,:],Xexp[trainset,:,:]], Yalt[trainset,:], batch_size=32)[1]
#Scores[1,0] = model.evaluate([Xchr[testset,:,:],Xexp[testset,:,:]], Yalt[testset,:], batch_size=32)[1]
#Scores[2,0] = model.evaluate([Xchr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
#pd.DataFrame(Scores, columns=['Accuracy']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomEndoAll/prediction_scores.csv')

