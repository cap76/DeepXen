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
from keras.models import load_model

np.random.seed(0)

Output2 = genfromtxt('BW_SomiteEndo_All/allpeaks_labelled_cum_final.se.bed',delimiter='\t',dtype=None)

#Load in the list of genes
#onlist = genfromtxt('ChIP_SomiteEndo_All/OnMemFDRtGenes.txt',delimiter='\t',dtype=None)
#downlist = genfromtxt('ChIP_SomiteEndo_All/ReprogrammedDowntGenes.txt',delimiter='\t',dtype=None)
#offlist = genfromtxt('ChIP_SomiteEndo_All/OffMemFDRtGenes.txt',delimiter='\t',dtype=None)
#uplist = genfromtxt('ChIP_SomiteEndo_All/ReprogrammedUptGenes.txt',delimiter='\t',dtype=None)
#complist = genfromtxt('ChIP_SomiteEndo_All/ComplementarySettGenes.txt',delimiter='\t',dtype=None)

#Generate the Y matrix
#Yalt = np.zeros((np.shape(Output2)[0],5))
#Xexpa = np.zeros((np.shape(Output2)[0],1,2))
#Xrawexp = np.zeros((np.shape(Output2)[0],7))

#for i in range(0, np.shape(Yalt)[0]):
#      Yalt[i,0] = (np.intersect1d(Output2['f3'][i],onlist)).size
#      Yalt[i,1] = (np.intersect1d(Output2['f3'][i],downlist)).size
#      Yalt[i,2] = (np.intersect1d(Output2['f3'][i],offlist)).size
#      Yalt[i,3] = (np.intersect1d(Output2['f3'][i],uplist)).size
#      Yalt[i,4] = 1-max(Yalt[i,0:4]) 
#
#explist = genfromtxt('ChIP_SomiteEndo_All/edger_de_pairing_BIGtable_endoIVF_somiteDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)
#explist1 = genfromtxt('ChIP_SomiteEndo_All/edger_de_pairing_BIGtable_endoIVF_somiteDonor_endoNT_plus_DE.nan.csv',delimiter='\t',dtype=None)
#
#for i in range(0, np.shape(Yalt)[0]):
#   arr_ind = np.where(Output2['f3'][i] == explist['f0'])
#   Xexpa[i,0,0] = np.log2(explist['f1'][arr_ind].astype('f')+1)
#   Xexpa[i,0,1] = np.log2(explist['f2'][arr_ind].astype('f')+1)

#np.nan_to_num(Xexpa,copy=False)
#
fil=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_SomiteEndo_All/*bwe.tab'))

#Load in all the expression data
encoder = LabelEncoder()
meh = encoder.fit(Output2['f0'])
#Xchr = np.zeros((np.shape(Output2)[0],600,np.shape(fil)[0]))
#for k in range(0, np.shape(fil)[0]):
#        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
#        Xchr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output

#Xexp = copy.deepcopy(Xexpa)
#Xexp[:,0,0] = ( Xexp[:,0,0] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0])
#Xexp[:,0,1] = ( Xexp[:,0,1] - np.nanmean(Xexp[trainset,0,1]) ) / np.nanstd(Xexp[trainset,0,1])

#Xchra = copy.deepcopy(Xchr)

#np.save('ResultsSomiteEndoAll/Xexp.npy', Xexp)
#np.save('ResultsSomiteEndoAll/Xexpa.npy', Xexpa)
#np.save('ResultsSomiteEndoAll/Xchr.npy', Xchr)

Xexp = np.load('ResultsSomiteEndoAll/Xexp.npy')
Xexpa = np.load('ResultsSomiteEndoAll/Xexpa.npy')
Xchr = np.load('ResultsSomiteEndoAll/Xchr.npy')

#Split into train/val/test sets
trainset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-18]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-17]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-16]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-15]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-14]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-13]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-12]))
valset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-11]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-10]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-9]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-8]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-7]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-6]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-5]))
testset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]))

#Xchr = np.load('ResultsSomiteEndoAll/Xchr.npy')

#Normalise histone modification data
for k in range(0, np.shape(fil)[0]):
      Xchr[:,0:600,k] = ( Xchr[:,0:600,k] - np.nanmean(Xchr[trainset,0:600,k]) ) / np.nanstd(Xchr[trainset,0:600,k])

plt.switch_backend('agg')
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
plt.switch_backend('agg')
from sklearn.decomposition import PCA

pca = PCA(n_components=20)
PCA1 = pca.fit_transform( np.concatenate( (np.mean(Xchr,axis=1),np.squeeze(Xexp)),axis=1))
np.save('ResultsSomiteEndoAll/PCA1.npy', PCA1)

PCA2 = pca.fit_transform( np.mean(Xchr,axis=1) )
np.save('ResultsSomiteEndoAll/PCA2.npy', PCA2)

import umap
reducer = umap.UMAP()
embedding1 = reducer.fit_transform(PCA1)
np.save('EndoFull/UMAP1.npy', embedding1)

embedding2 = reducer.fit_transform(  )
np.save('EndoFull/UMAP2.npy', embedding2)



