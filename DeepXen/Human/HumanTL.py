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

np.random.seed(0)

Output2 = genfromtxt('BWCM/EndoNT_On_all.inde.bed',delimiter='\t',dtype=None)
#Load in the list of genes
onlist = genfromtxt('CM/ON_HDF_IVF_NT_CM.csv',delimiter='\t',dtype=None)
downlist = genfromtxt('CM/RD_HDF_IVF_NT_CM.csv',delimiter='\t',dtype=None)
offlist = genfromtxt('CM/OFF_HDF_IVF_NT_CM.csv',delimiter='\t',dtype=None)
uplist = genfromtxt('CM/RU_HDF_IVF_NT_CM.csv',delimiter='\t',dtype=None)
complist = genfromtxt('CM/ComplementarySetGenes.txt',delimiter='\t',dtype=None)
#Generate the Y matrix
Yalt = np.zeros((np.shape(Output2)[0],5))
Xexpa = np.zeros((np.shape(Output2)[0],1,2))
for i in range(0, np.shape(Yalt)[0]):
      Yalt[i,0] = (np.intersect1d(Output2['f3'][i],onlist)).size
      Yalt[i,1] = (np.intersect1d(Output2['f3'][i],downlist)).size
      Yalt[i,2] = (np.intersect1d(Output2['f3'][i],offlist)).size
      Yalt[i,3] = (np.intersect1d(Output2['f3'][i],uplist)).size
      Yalt[i,4] = 1-max(Yalt[i,0:4]) 

explist = genfromtxt('FPKM.csv',delimiter='\t',dtype=None, skip_header=1)

for i in range(0, np.shape(Yalt)[0]):
   arr_ind = np.where(Output2['f3'][i] == explist['f0'])
   x1 = (explist['f1'][arr_ind].astype('f')+explist['f2'][arr_ind].astype('f')+explist['f3'][arr_ind].astype('f'))/3
   x2 = (explist['f4'][arr_ind].astype('f')+explist['f5'][arr_ind].astype('f')+explist['f10'][arr_ind].astype('f')+explist['f11'][arr_ind].astype('f'))/4
   Xexpa[i,0,0] = np.log2(x1+1)
   #explist['f1'][arr_ind].astype('f')+1)
   Xexpa[i,0,1] = np.log2(x2+1)
   #explist['f2'][arr_ind].astype('f')+1)

np.nan_to_num(Xexpa,copy=False)
fil=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/Human/BWCM/*.tab'))
#Load in all the expression data
encoder = LabelEncoder()
meh = encoder.fit(Output2['f0'])
Xchr = np.zeros((np.shape(Output2)[0],600,np.shape(fil)[0]))
for k in range(0, np.shape(fil)[0]):
        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
        Xchr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output

#Chromatin mods
for k in range(0, np.shape(fil)[0]):
        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
        Xchr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output

#Split into train/val/test sets
#trainset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-18]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-17]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-16]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-15]) | (Output2['f0']==meh.classes_[np.shape(meh.classes$
#valset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-11]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-10]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-9]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-8]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0$
#testset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]))
trainset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-18]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-17]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-16]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-15]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-14]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-13]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-12]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-22]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-23]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-24]))
valset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-11]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-10]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-9]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-8]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-7]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-6]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-5]))
testset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-19]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-20]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-21]))



#Normalise the expression data
Xexp = copy.deepcopy(Xexpa)
Xexp[:,0,0] = ( Xexp[:,0,0] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0])
Xexp[:,0,1] = ( Xexp[:,0,1] - np.nanmean(Xexp[trainset,0,1]) ) / np.nanstd(Xexp[trainset,0,1])

Xchra = copy.deepcopy(Xchr)

#Normalise histone modification data
for k in range(0, np.shape(fil)[0]):
      Xchr[:,0:600,k] = ( Xchr[:,0:600,k] - np.nanmean(Xchr[trainset,0:600,k]) ) / np.nanstd(Xchr[trainset,0:600,k])


#Now we can load in the other two datasets and for predictions. 

#Finally load in the model ...
Output3 = genfromtxt('BWEC/EndoNT_On_all.inde.bed',delimiter='\t',dtype=None)

#Load in the list of genes
onlist3 = genfromtxt('EC/ON_HDF_IVF_NT_EC.csv',delimiter='\t',dtype=None)
downlist3 = genfromtxt('EC/RD_HDF_IVF_NT_EC.csv',delimiter='\t',dtype=None)
offlist3 = genfromtxt('EC/OFF_HDF_IVF_NT_EC.csv',delimiter='\t',dtype=None)
uplist3 = genfromtxt('EC/RU_HDF_IVF_NT_EC.csv',delimiter='\t',dtype=None)
complist3 = genfromtxt('EC/ComplementarySetGenes.txt',delimiter='\t',dtype=None)

encoder = LabelEncoder()
meh = encoder.fit(Output3['f0'])
trainset3=((Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-18]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-17]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-16]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-15]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-14]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-13]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-12]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-22]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-23]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-24]))
valset3=((Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-11]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-10]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-9]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-8]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-7]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-6]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-5]))
testset3=((Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-19]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-20]) | (Output3['f0']==meh.classes_[np.shape(meh.classes_)[0]-21]))


#Generate the Y matrix
Yalt3 = np.zeros((np.shape(Output3)[0],5))
Xexpa3 = np.zeros((np.shape(Output3)[0],1,2))

for i in range(0, np.shape(Yalt3)[0]):
      Yalt3[i,0] = (np.intersect1d(Output3['f3'][i],onlist3)).size
      Yalt3[i,1] = (np.intersect1d(Output3['f3'][i],downlist3)).size
      Yalt3[i,2] = (np.intersect1d(Output3['f3'][i],offlist3)).size
      Yalt3[i,3] = (np.intersect1d(Output3['f3'][i],uplist3)).size
      Yalt3[i,4] = 1-max(Yalt3[i,0:4]) 

for i in range(0, np.shape(Yalt3)[0]):
   arr_ind = np.where(Output3['f3'][i] == explist['f0'])
   x1 = (explist['f1'][arr_ind].astype('f')+explist['f2'][arr_ind].astype('f')+explist['f3'][arr_ind].astype('f'))/3
   x2 = (explist['f6'][arr_ind].astype('f')+explist['f7'][arr_ind].astype('f')+explist['f12'][arr_ind].astype('f')+explist['f13'][arr_ind].astype('f'))/4
   Xexpa3[i,0,0] = np.log2(x1+1)
   Xexpa3[i,0,1] = np.log2(x2+1)

np.nan_to_num(Xexpa3,copy=False)

fil=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/Human/BWEC/*.tab'))

#Load in all the expression data
encoder = LabelEncoder()
meh = encoder.fit(Output2['f0'])
Xchr3 = np.zeros((np.shape(Output3)[0],600,np.shape(fil)[0]))
for k in range(0, np.shape(fil)[0]):
        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
        Xchr3[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output

#Normalise the expression data based on training dataset only
Xexp3 = copy.deepcopy(Xexpa3)
Xexp3[:,0,0] = ( Xexp3[:,0,0] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0])
Xexp3[:,0,1] = ( Xexp3[:,0,1] - np.nanmean(Xexp[trainset,0,1]) ) / np.nanstd(Xexp[trainset,0,1])

Xchra3 = copy.deepcopy(Xchr3)

#Normalise histone modification data
for k in range(0, np.shape(fil)[0]):
      Xchr3[:,0:600,k] = ( Xchr3[:,0:600,k] - np.nanmean(Xchr[trainset,0:600,k]) ) / np.nanstd(Xchr[trainset,0:600,k])




#Load ESC dataset

Output4 = genfromtxt('BWESC/EndoNT_On_all.inde.bed',delimiter='\t',dtype=None)

#Load in the list of genes
onlist4 = genfromtxt('ESC/ON_HDF_IVF_NT_ESC.csv',delimiter='\t',dtype=None)
downlist4 = genfromtxt('ESC/RD_HDF_IVF_NT_ESC.csv',delimiter='\t',dtype=None)
offlist4 = genfromtxt('ESC/OFF_HDF_IVF_NT_ESC.csv',delimiter='\t',dtype=None)
uplist4 = genfromtxt('ESC/RU_HDF_IVF_NT_ESC.csv',delimiter='\t',dtype=None)
complist4 = genfromtxt('ESC/ComplementarySetGenes.txt',delimiter='\t',dtype=None)

encoder = LabelEncoder()
meh = encoder.fit(Output4['f0'])
trainset4=((Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-18]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-17]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-16]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-15]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-14]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-13]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-12]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-22]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-23]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-24]))
valset4=((Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-11]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-10]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-9]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-8]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-7]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-6]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-5]))
testset4=((Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-19]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-20]) | (Output4['f0']==meh.classes_[np.shape(meh.classes_)[0]-21]))



#Generate the Y matrix
Yalt4 = np.zeros((np.shape(Output4)[0],5))
Xexpa4 = np.zeros((np.shape(Output4)[0],1,2))

for i in range(0, np.shape(Yalt4)[0]):
      Yalt4[i,0] = (np.intersect1d(Output4['f3'][i],onlist)).size
      Yalt4[i,1] = (np.intersect1d(Output4['f3'][i],downlist)).size
      Yalt4[i,2] = (np.intersect1d(Output4['f3'][i],offlist)).size
      Yalt4[i,3] = (np.intersect1d(Output4['f3'][i],uplist)).size
      Yalt4[i,4] = 1-max(Yalt4[i,0:4]) 

for i in range(0, np.shape(Yalt4)[0]):
   arr_ind = np.where(Output4['f3'][i] == explist['f0'])
   x1 = (explist['f1'][arr_ind].astype('f')+explist['f2'][arr_ind].astype('f')+explist['f3'][arr_ind].astype('f'))/3
   x2 = (explist['f8'][arr_ind].astype('f')+explist['f9'][arr_ind].astype('f')+explist['f14'][arr_ind].astype('f')+explist['f15'][arr_ind].astype('f'))/4
   Xexpa4[i,0,0] = np.log2(x1+1)
   Xexpa4[i,0,1] = np.log2(x2+1)

np.nan_to_num(Xexpa4,copy=False)

fil=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/Human/BWESC/*.tab'))

#Load in all the expression data
encoder = LabelEncoder()
meh = encoder.fit(Output4['f0'])
Xchr4 = np.zeros((np.shape(Output4)[0],600,np.shape(fil)[0]))
for k in range(0, np.shape(fil)[0]):
        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
        Xchr4[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output

#Normalise the expression data
Xexp4 = copy.deepcopy(Xexpa4)
Xexp4[:,0,0] = ( Xexp4[:,0,0] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0])
Xexp4[:,0,1] = ( Xexp4[:,0,1] - np.nanmean(Xexp[trainset,0,1]) ) / np.nanstd(Xexp[trainset,0,1])

Xchra4 = copy.deepcopy(Xchr4)

#Normalise histone modification data
for k in range(0, np.shape(fil)[0]):
      Xchr4[:,0:600,k] = ( Xchr4[:,0:600,k] - np.nanmean(Xchr[trainset,0:600,k]) ) / np.nanstd(Xchr[trainset,0:600,k])

#Now load in the models ...
from keras.models import load_model
model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/Human/BWCM/Model_deeper.h5') 

from sklearn.metrics import roc_curve, precision_recall_curve,auc

fpr_keras, tpr_keras, thresholds_keras = roc_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())

fpr_keras3, tpr_keras3, thresholds_keras = roc_curve(Yalt3.ravel(), model.predict([Xchr3,Xexp3], batch_size=1000).ravel())
pr_keras3, re_keras3, thresholds_keras = precision_recall_curve(Yalt3.ravel(), model.predict([Xchr3,Xexp3], batch_size=1000).ravel())

fpr_keras4, tpr_keras4, thresholds_keras = roc_curve(Yalt4.ravel(), model.predict([Xchr4,Xexp4], batch_size=1000).ravel())
pr_keras4, re_keras4, thresholds_keras = precision_recall_curve(Yalt4.ravel(), model.predict([Xchr4,Xexp4], batch_size=1000).ravel())

np.save('TL/fpr_CM_', fpr_keras)
np.save('TL/tpr_CM_', tpr_keras)
np.save('TL/pr_CM_', pr_keras)
np.save('TL/re_CM_', re_keras)

np.save('TL/fpr_EC_', fpr_keras3)
np.save('TL/tpr_EC_', tpr_keras3)
np.save('TL/pr_EC_', pr_keras3)
np.save('TL/re_EC_', re_keras3)

np.save('TL/fpr_CM_ESC_', fpr_keras4)
np.save('TL/tpr_CM_ESC_', tpr_keras4)
np.save('TL/pr_CM_ESC_', pr_keras4)
np.save('TL/re_CM_ESC_', re_keras4)


for layer in model.layers[:25]:
   layer.trainable = False

nclasses = np.sum(Yalt3[trainset3,:],0)
cwe = np.ndarray.max(nclasses) / nclasses
class_weight={0: cwe[0], 1: cwe[1], 2: cwe[2], 3: cwe[3], 4: cwe[4]} 
sgd = optimizers.SGD(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss="categorical_crossentropy", optimizer = sgd, metrics=['categorical_accuracy'])
callbacks = [GetBest(monitor='val_categorical_accuracy', verbose=1, mode='max')]
model.fit([Xchr3[trainset3,:,:],Xexp3[trainset3,:,:]], Yalt3[trainset3,:], validation_data=([Xchr3[testset3,:,:],Xexp3[testset3,:,:]], Yalt3[testset3,:]), epochs=500, batch_size=1000,class_weight = class_weight,callbacks=callbacks)
model.save('TL/CM_to_EC_TL.h5')

fpr_keras5, tpr_keras5, thresholds_keras5 = roc_curve(Yalt3[valset3,:].ravel(), model.predict([Xchr3[valset3,:,:],Xexp3[valset3,:,:]], batch_size=1000).ravel())
pr_keras5, re_keras5, thresholds_keras5 = precision_recall_curve(Yalt3[valset3,:].ravel(), model.predict([Xchr3[valset3,:,:],Xexp3[valset3,:,:]], batch_size=1000).ravel())

np.save('TL/fpr_CM_EC_TL', fpr_keras5)
np.save('TL/tpr_CM_EC_TL', tpr_keras5)
np.save('TL/pr_CM_EC_TL', pr_keras5)
np.save('TL/re_CM_EC_TL', re_keras5)

from keras.models import load_model
model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/Human/BWCM/Model_deeper.h5') 

for layer in model.layers[:25]:
   layer.trainable = False

nclasses = np.sum(Yalt4[trainset4,:],0)
cwe = np.ndarray.max(nclasses) / nclasses
class_weight={0: cwe[0], 1: cwe[1], 2: cwe[2], 3: cwe[3], 4: cwe[4]} 
sgd = optimizers.SGD(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss="categorical_crossentropy", optimizer = sgd, metrics=['categorical_accuracy'])
callbacks = [GetBest(monitor='val_categorical_accuracy', verbose=1, mode='max')]
model.fit([Xchr4[trainset4,:,:],Xexp4[trainset4,:,:]], Yalt4[trainset4,:], validation_data=([Xchr4[testset4,:,:],Xexp4[testset4,:,:]], Yalt4[testset4,:]), epochs=500, batch_size=1000,class_weight = class_weight,callbacks=callbacks)
model.save('TL/CM_to_ESC_TL.h5')

fpr_keras6, tpr_keras6, thresholds_keras5 = roc_curve(Yalt4[valset4,:].ravel(), model.predict([Xchr4[valset4,:,:],Xexp4[valset4,:,:]], batch_size=1000).ravel())
pr_keras6, re_keras6, thresholds_keras5 = precision_recall_curve(Yalt4[valset4,:].ravel(), model.predict([Xchr4[valset4,:,:],Xexp4[valset4,:,:]], batch_size=1000).ravel())

np.save('TL/fpr_CM_EC_TL', fpr_keras6)
np.save('TL/tpr_CM_EC_TL', tpr_keras6)
np.save('TL/pr_CM_EC_TL', pr_keras6)
np.save('TL/re_CM_EC_TL', re_keras6)
