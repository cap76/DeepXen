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

fil1=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/Endo*bwe.tab')) 
fil2=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/Ecto*bwe.tab'))
fil3=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/GC.bwe.tab'))
fil4=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/*Methylation.bwe.tab'))

fil=fil1+fil2+fil3+fil4
filr=fil2+fil1+fil3+fil4


#Load in all the expression data
encoder = LabelEncoder()
meh = encoder.fit(Output2['f0'])
Xchr = np.zeros((np.shape(Output2)[0],600,np.shape(fil)[0]))
Xchrr = np.zeros((np.shape(Output2)[0],600,np.shape(fil)[0]))
for k in range(0, np.shape(fil)[0]):
        Output = genfromtxt(fil[k],delimiter='\t',dtype=None, skip_header=3)
        Xchr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output
        Output = genfromtxt(filr[k],delimiter='\t',dtype=None, skip_header=3)
        Xchrr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output



#Split into train/val/test sets
trainset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-18]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-17]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-16]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-15]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-14]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-13]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-12]))
valset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-11]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-10]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-9]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-8]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-7]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-6]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-5]))
testset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]))

#Normalise the expression data
Xexp = copy.deepcopy(Xexpa)
Xexp[:,0,0] = ( Xexp[:,0,0] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0])
Xexp[:,0,1] = ( Xexp[:,0,1] - np.nanmean(Xexp[trainset,0,1]) ) / np.nanstd(Xexp[trainset,0,1])

Xchra = copy.deepcopy(Xchr)

#Normalise histone modification data
for k in range(0, np.shape(fil)[0]):
      Xchr[:,0:600,k] = ( Xchr[:,0:600,k] - np.nanmean(Xchr[trainset,0:600,k]) ) / np.nanstd(Xchr[trainset,0:600,k])
      Xchrr[:,0:600,k] = ( Xchrr[:,0:600,k] - np.nanmean(Xchrr[trainset,0:600,k]) ) / np.nanstd(Xchrr[trainset,0:600,k])

from keras.models import load_model
model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAll/Model_deeper.h5')
model1 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsSomiteEndoAll/Model_deeper_seed=2.h5')
model2 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_TL.h5')
model3 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_TL_seed=2.h5')

#predictions1 = model.predict([Xchr,Xexp], batch_size=1000)
#predictions2 = model1.predict([Xchr,Xexp], batch_size=1000)
#predictions3 = model2.predict([Xchr,Xexp], batch_size=1000)
#predictions4 = model3.predict([Xchr,Xexp], batch_size=1000)

#plt.switch_backend('agg')

Scores = np.zeros((6,1))
Scores[0,0] = model.evaluate([Xchr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
Scores[1,0] = model1.evaluate([Xchr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
Scores[2,0] = model.evaluate([Xchrr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
Scores[3,0] = model1.evaluate([Xchrr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]

Scores[4,0] = model2.evaluate([Xchr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
Scores[5,0] = model3.evaluate([Xchr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]

pd.DataFrame(Scores, columns=['Accuracy']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/PredictionScores_TL.csv')

#df1 = pd.DataFrame(Output2['f0'],columns=['chr'])
#df2 = pd.DataFrame(Output2['f1'],columns=['start'])
#df3 = pd.DataFrame(Output2['f2'],columns=['end'])
#df4 = pd.DataFrame(Output2['f3'],columns=['C1'])
#df5 = pd.DataFrame(Output2['f4'],columns=['C2'])
#df6 = pd.DataFrame(Output2['f5'],columns=['C3'])
#df7 = pd.DataFrame(Output2['f6'],columns=['C3'])
#df8 = pd.DataFrame(Output2['f7'],columns=['C4'])
#df9 = pd.DataFrame(Output2['f8'],columns=['C5'])
#df10 = pd.DataFrame(Output2['f9'],columns=['C6'])
#df11 = pd.DataFrame(Output2['f10'],columns=['C7'])
#df12 = pd.DataFrame(Yalt[:,0],columns=['C8'])
#df13 = pd.DataFrame(Yalt[:,1],columns=['C9'])
#df14 = pd.DataFrame(Yalt[:,2],columns=['C10'])
#df15 = pd.DataFrame(Yalt[:,3],columns=['C11'])
#df16 = pd.DataFrame(Yalt[:,4],columns=['C12'])

#df14b = pd.DataFrame(predictions,columns=['pC1','pC2','pC3','pC4','pC5'])

#prediction = np.concatenate((df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14,df15,df16, df14b),1)
#pd.DataFrame(prediction, columns=['Chr','start','end','Gene','IVF','Donor','NT','FC','pVal','FC','pVal','ON','RepDown','Off','RepUp','Neg','pC1','pC2','pC3','pC4','pC5']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/prediction_deeper_TL.csv')
#NOW GENERATE SOME EXAMPLE PLOTS
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
plt.switch_backend('agg')

from sklearn.metrics import roc_curve, precision_recall_curve,auc

fpr_keras, tpr_keras, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())

fpr_keras1, tpr_keras1, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model1.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_keras1, re_keras1, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model1.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())

fpr_keras2, tpr_keras2, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model.predict([Xchrr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_keras2, re_keras2, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model.predict([Xchrr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())

fpr_keras3, tpr_keras3, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model1.predict([Xchrr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_keras3, re_keras3, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model1.predict([Xchrr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())


fpr_keras5, tpr_keras5, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model2.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_keras5, re_keras5, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model2.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())

fpr_keras6, tpr_keras6, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model3.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_keras6, re_keras6, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model3.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())

#fpr_keras6, tpr_keras6, thresholds_keras = roc_curve(Yalt[valset,:].ravel(), model3.predict(Xchr[valset,:,:], batch_size=1000).ravel())
#pr_keras6, re_keras6, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model3.predict(Xchr[valset,:,:], batch_size=1000).ravel())

Y0 = Yalt[valset,:]
P0 = model.predict([Xchr[valset,:,:],Xexp[valset,:,:]])
Y1 = Yalt[valset,:]
P1 = model1.predict([Xchr[valset,:,:],Xexp[valset,:,:]])
Y2 = Yalt[valset,:]
P2 = model2.predict([Xchr[valset,:,:],Xexp[valset,:,:]])
Y3 = Yalt[valset,:]
P3 = model3.predict([Xchr[valset,:,:],Xexp[valset,:,:]])

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


fpr_keras2_1, tpr_keras2_1, thresholds_kerast = roc_curve(Y2[:,[True,False,False,False,False]], P2[:,[True,False,False,False,False]])
fpr_keras2_2, tpr_keras2_2, thresholds_kerast = roc_curve(Y2[:,[False,True,False,False,False]], P2[:,[False,True,False,False,False]])
fpr_keras2_3, tpr_keras2_3, thresholds_kerast = roc_curve(Y2[:,[False,False,True,False,False]], P2[:,[False,False,True,False,False]])
fpr_keras2_4, tpr_keras2_4, thresholds_kerast = roc_curve(Y2[:,[False,False,False,True,False]], P2[:,[False,False,False,True,False]])
fpr_keras2_5, tpr_keras2_5, thresholds_kerast = roc_curve(Y2[:,[False,False,False,False,True]], P2[:,[False,False,False,False,True]])

fpr_keras3_1, tpr_keras3_1, thresholds_kerast = roc_curve(Y3[:,[True,False,False,False,False]], P3[:,[True,False,False,False,False]])
fpr_keras3_2, tpr_keras3_2, thresholds_kerast = roc_curve(Y3[:,[False,True,False,False,False]], P3[:,[False,True,False,False,False]])
fpr_keras3_3, tpr_keras3_3, thresholds_kerast = roc_curve(Y3[:,[False,False,True,False,False]], P3[:,[False,False,True,False,False]])
fpr_keras3_4, tpr_keras3_4, thresholds_kerast = roc_curve(Y3[:,[False,False,False,True,False]], P3[:,[False,False,False,True,False]])
fpr_keras3_5, tpr_keras3_5, thresholds_kerast = roc_curve(Y3[:,[False,False,False,False,True]], P3[:,[False,False,False,False,True]])

np.save('ResultsEndoAll/tpr_kears0_1_endoTargOnly_TL.npy', tpr_keras0_1)
np.save('ResultsEndoAll/tpr_kears0_2_endoTargOnly_TL.npy', tpr_keras0_2)
np.save('ResultsEndoAll/tpr_kears0_3_endoTargOnly_TL.npy', tpr_keras0_3)
np.save('ResultsEndoAll/tpr_kears0_4_endoTargOnly_TL.npy', tpr_keras0_4)
np.save('ResultsEndoAll/tpr_kears0_5_endoTargOnly_TL.npy', tpr_keras0_5)

np.save('ResultsEndoAll/fpr_kears0_1_endoTargOnly_TL.npy', fpr_keras0_1)
np.save('ResultsEndoAll/fpr_kears0_2_endoTargOnly_TL.npy', fpr_keras0_2)
np.save('ResultsEndoAll/fpr_kears0_3_endoTargOnly_TL.npy', fpr_keras0_3)
np.save('ResultsEndoAll/fpr_kears0_4_endoTargOnly_TL.npy', fpr_keras0_4)
np.save('ResultsEndoAll/fpr_kears0_5_endoTargOnly_TL.npy', fpr_keras0_5)

np.save('ResultsEndoAll/tpr_kears1_1_endoTargOnly_TL.npy', tpr_keras1_1)
np.save('ResultsEndoAll/tpr_kears1_2_endoTargOnly_TL.npy', tpr_keras1_2)
np.save('ResultsEndoAll/tpr_kears1_3_endoTargOnly_TL.npy', tpr_keras1_3)
np.save('ResultsEndoAll/tpr_kears1_4_endoTargOnly_TL.npy', tpr_keras1_4)
np.save('ResultsEndoAll/tpr_kears1_5_endoTargOnly_TL.npy', tpr_keras1_5)

np.save('ResultsEndoAll/fpr_kears1_1_endoTargOnly_TL.npy', fpr_keras1_1)
np.save('ResultsEndoAll/fpr_kears1_2_endoTargOnly_TL.npy', fpr_keras1_2)
np.save('ResultsEndoAll/fpr_kears1_3_endoTargOnly_TL.npy', fpr_keras1_3)
np.save('ResultsEndoAll/fpr_kears1_4_endoTargOnly_TL.npy', fpr_keras1_4)
np.save('ResultsEndoAll/fpr_kears1_5_endoTargOnly_TL.npy', fpr_keras1_5)

np.save('ResultsEndoAll/tpr_kears2_1_endoTargOnly_TL.npy', tpr_keras2_1)
np.save('ResultsEndoAll/tpr_kears2_2_endoTargOnly_TL.npy', tpr_keras2_2)
np.save('ResultsEndoAll/tpr_kears2_3_endoTargOnly_TL.npy', tpr_keras2_3)
np.save('ResultsEndoAll/tpr_kears2_4_endoTargOnly_TL.npy', tpr_keras2_4)
np.save('ResultsEndoAll/tpr_kears2_5_endoTargOnly_TL.npy', tpr_keras2_5)

np.save('ResultsEndoAll/fpr_kears2_1_endoTargOnly_TL.npy', fpr_keras2_1)
np.save('ResultsEndoAll/fpr_kears2_2_endoTargOnly_TL.npy', fpr_keras2_2)
np.save('ResultsEndoAll/fpr_kears2_3_endoTargOnly_TL.npy', fpr_keras2_3)
np.save('ResultsEndoAll/fpr_kears2_4_endoTargOnly_TL.npy', fpr_keras2_4)
np.save('ResultsEndoAll/fpr_kears2_5_endoTargOnly_TL.npy', fpr_keras2_5)

np.save('ResultsEndoAll/tpr_kears3_1_endoTargOnly_TL.npy', tpr_keras3_1)
np.save('ResultsEndoAll/tpr_kears3_2_endoTargOnly_TL.npy', tpr_keras3_2)
np.save('ResultsEndoAll/tpr_kears3_3_endoTargOnly_TL.npy', tpr_keras3_3)
np.save('ResultsEndoAll/tpr_kears3_4_endoTargOnly_TL.npy', tpr_keras3_4)
np.save('ResultsEndoAll/tpr_kears3_5_endoTargOnly_TL.npy', tpr_keras3_5)

np.save('ResultsEndoAll/fpr_kears3_1_endoTargOnly_TL.npy', fpr_keras3_1)
np.save('ResultsEndoAll/fpr_kears3_2_endoTargOnly_TL.npy', fpr_keras3_2)
np.save('ResultsEndoAll/fpr_kears3_3_endoTargOnly_TL.npy', fpr_keras3_3)
np.save('ResultsEndoAll/fpr_kears3_4_endoTargOnly_TL.npy', fpr_keras3_4)
np.save('ResultsEndoAll/fpr_kears3_5_endoTargOnly_TL.npy', fpr_keras3_5)




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

pr_keras2_1, re_keras2_1, thresholds_kerast = precision_recall_curve(Y2[:,[True,False,False,False,False]], P2[:,[True,False,False,False,False]])
pr_keras2_2, re_keras2_2, thresholds_kerast = precision_recall_curve(Y2[:,[False,True,False,False,False]], P2[:,[False,True,False,False,False]])
pr_keras2_3, re_keras2_3, thresholds_kerast = precision_recall_curve(Y2[:,[False,False,True,False,False]], P2[:,[False,False,True,False,False]])
pr_keras2_4, re_keras2_4, thresholds_kerast = precision_recall_curve(Y2[:,[False,False,False,True,False]], P2[:,[False,False,False,True,False]])
pr_keras2_5, re_keras2_5, thresholds_kerast = precision_recall_curve(Y2[:,[False,False,False,False,True]], P2[:,[False,False,False,False,True]])

pr_keras3_1, re_keras3_1, thresholds_kerast = precision_recall_curve(Y3[:,[True,False,False,False,False]], P3[:,[True,False,False,False,False]])
pr_keras3_2, re_keras3_2, thresholds_kerast = precision_recall_curve(Y3[:,[False,True,False,False,False]], P3[:,[False,True,False,False,False]])
pr_keras3_3, re_keras3_3, thresholds_kerast = precision_recall_curve(Y3[:,[False,False,True,False,False]], P3[:,[False,False,True,False,False]])
pr_keras3_4, re_keras3_4, thresholds_kerast = precision_recall_curve(Y3[:,[False,False,False,True,False]], P3[:,[False,False,False,True,False]])
pr_keras3_5, re_keras3_5, thresholds_kerast = precision_recall_curve(Y3[:,[False,False,False,False,True]], P3[:,[False,False,False,False,True]])


np.save('ResultsEndoAll/pr_kears0_1_endoTargOnly_TL.npy', pr_keras0_1)
np.save('ResultsEndoAll/pr_kears0_2_endoTargOnly_TL.npy', pr_keras0_2)
np.save('ResultsEndoAll/pr_kears0_3_endoTargOnly_TL.npy', pr_keras0_3)
np.save('ResultsEndoAll/pr_kears0_4_endoTargOnly_TL.npy', pr_keras0_4)
np.save('ResultsEndoAll/pr_kears0_5_endoTargOnly_TL.npy', pr_keras0_5)

np.save('ResultsEndoAll/re_kears0_1_endoTargOnly_TL.npy', re_keras0_1)
np.save('ResultsEndoAll/re_kears0_2_endoTargOnly_TL.npy', re_keras0_2)
np.save('ResultsEndoAll/re_kears0_3_endoTargOnly_TL.npy', re_keras0_3)
np.save('ResultsEndoAll/re_kears0_4_endoTargOnly_TL.npy', re_keras0_4)
np.save('ResultsEndoAll/re_kears0_5_endoTargOnly_TL.npy', re_keras0_5)


np.save('ResultsEndoAll/pr_kears1_1_endoTargOnly_TL.npy', pr_keras1_1)
np.save('ResultsEndoAll/pr_kears1_2_endoTargOnly_TL.npy', pr_keras1_2)
np.save('ResultsEndoAll/pr_kears1_3_endoTargOnly_TL.npy', pr_keras1_3)
np.save('ResultsEndoAll/pr_kears1_4_endoTargOnly_TL.npy', pr_keras1_4)
np.save('ResultsEndoAll/pr_kears1_5_endoTargOnly_TL.npy', pr_keras1_5)

np.save('ResultsEndoAll/re_kears1_1_endoTargOnly_TL.npy', re_keras1_1)
np.save('ResultsEndoAll/re_kears1_2_endoTargOnly_TL.npy', re_keras1_2)
np.save('ResultsEndoAll/re_kears1_3_endoTargOnly_TL.npy', re_keras1_3)
np.save('ResultsEndoAll/re_kears1_4_endoTargOnly_TL.npy', re_keras1_4)
np.save('ResultsEndoAll/re_kears1_5_endoTargOnly_TL.npy', re_keras1_5)

np.save('ResultsEndoAll/pr_kears2_1_endoTargOnly_TL.npy', pr_keras2_1)
np.save('ResultsEndoAll/pr_kears2_2_endoTargOnly_TL.npy', pr_keras2_2)
np.save('ResultsEndoAll/pr_kears2_3_endoTargOnly_TL.npy', pr_keras2_3)
np.save('ResultsEndoAll/pr_kears2_4_endoTargOnly_TL.npy', pr_keras2_4)
np.save('ResultsEndoAll/pr_kears2_5_endoTargOnly_TL.npy', pr_keras2_5)

np.save('ResultsEndoAll/re_kears2_1_endoTargOnly_TL.npy', re_keras2_1)
np.save('ResultsEndoAll/re_kears2_2_endoTargOnly_TL.npy', re_keras2_2)
np.save('ResultsEndoAll/re_kears2_3_endoTargOnly_TL.npy', re_keras2_3)
np.save('ResultsEndoAll/re_kears2_4_endoTargOnly_TL.npy', re_keras2_4)
np.save('ResultsEndoAll/re_kears2_5_endoTargOnly_TL.npy', re_keras2_5)


np.save('ResultsEndoAll/pr_kears3_1_endoTargOnly_TL.npy', pr_keras3_1)
np.save('ResultsEndoAll/pr_kears3_2_endoTargOnly_TL.npy', pr_keras3_2)
np.save('ResultsEndoAll/pr_kears3_3_endoTargOnly_TL.npy', pr_keras3_3)
np.save('ResultsEndoAll/pr_kears3_4_endoTargOnly_TL.npy', pr_keras3_4)
np.save('ResultsEndoAll/pr_kears3_5_endoTargOnly_TL.npy', pr_keras3_5)

np.save('ResultsEndoAll/re_kears3_1_endoTargOnly_TL.npy', re_keras3_1)
np.save('ResultsEndoAll/re_kears3_2_endoTargOnly_TL.npy', re_keras3_2)
np.save('ResultsEndoAll/re_kears3_3_endoTargOnly_TL.npy', re_keras3_3)
np.save('ResultsEndoAll/re_kears3_4_endoTargOnly_TL.npy', re_keras3_4)
np.save('ResultsEndoAll/re_kears3_5_endoTargOnly_TL.npy', re_keras3_5)




np.save('EndoFull/fpr_kears_TL.npy', fpr_keras)
np.save('EndoFull/fpr_kears1_TL.npy', fpr_keras1)
np.save('EndoFull/fpr_kears2_TL.npy', fpr_keras2)
np.save('EndoFull/fpr_kears3_TL.npy', fpr_keras3)
np.save('EndoFull/fpr_kears5_TL.npy', fpr_keras5)
np.save('EndoFull/fpr_kears6_TL.npy', fpr_keras6)

np.save('EndoFull/tpr_kears_TL.npy', tpr_keras)
np.save('EndoFull/tpr_kears1_TL.npy', tpr_keras1)
np.save('EndoFull/tpr_kears2_TL.npy', tpr_keras2)
np.save('EndoFull/tpr_kears3_TL.npy', tpr_keras3)
np.save('EndoFull/tpr_kears5_TL.npy', tpr_keras5)
np.save('EndoFull/tpr_kears6_TL.npy', tpr_keras6)

np.save('EndoFull/pr_kears_TL.npy', pr_keras)
np.save('EndoFull/pr_kears1_TL.npy', pr_keras1)
np.save('EndoFull/pr_kears2_TL.npy', pr_keras2)
np.save('EndoFull/pr_kears3_TL.npy', pr_keras3)
np.save('EndoFull/pr_kears5_TL.npy', pr_keras5)
np.save('EndoFull/pr_kears6_TL.npy', pr_keras6)

np.save('EndoFull/re_kears_TL.npy', re_keras)
np.save('EndoFull/re_kears1_TL.npy', re_keras1)
np.save('EndoFull/re_kears2_TL.npy', re_keras2)
np.save('EndoFull/re_kears3_TL.npy', re_keras3)
np.save('EndoFull/re_kears5_TL.npy', re_keras5)
np.save('EndoFull/re_kears6_TL.npy', re_keras6)


plt.figure()
plt.plot(re_keras5, pr_keras5, label='SomEndo') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(re_keras6, pr_keras6, label='SomEndo (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(re_keras2, pr_keras2, label='TL') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(re_keras3, pr_keras3, label='TL (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
#plt.xlim([0.0, 1.0])
#plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig('EndoFull/PR_TL.pdf')
plt.draw()

plt.figure()
plt.plot(fpr_keras5, tpr_keras5, label='SomEndo') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(fpr_keras6, tpr_keras6, label='SomEndo (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(fpr_keras2, tpr_keras2, label='TL') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(fpr_keras3, tpr_keras3, label='TL') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot([0, 1], [0, 1], 'k--', label='Random')
plt.plot([0,0, 1], [0, 1, 1], 'k-')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig('EndoFull/AUC_TL.pdf')
plt.draw()



classpred1 = model2.predict([Xchr,Xexp], batch_size=1000)
classpred2 = model3.predict([Xchr,Xexp], batch_size=1000)

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

M = np.zeros((6,26))
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

classpred = classpred2
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


M[1,0]= sum(OnP) #,val=True)
M[1,1]= sum(OnT) #,val=True) #.shape[0]
M[1,2]= sum(OnTP) #,val=True) #.shape[0]
M[1,3]= sum(OnFP) #,val=True) #.shape[0]
M[1,4]= sum(OnTN) #,val=True) #.shape[0]
M[1,5]= sum(OnFN) #,val=True) #.shape[0]
M[1,6]= sum(OffP) #,val=True) #.shape[0]
M[1,7]= sum(OffT) #,val=True) #.shape[0]
M[1,8]= sum(OffTP) #,val=True) #.shape[0]
M[1,9]= sum(OffFP) #,val=True) #.shape[0]
M[1,10]= sum(OnTN) #,val=True) #.shape[0]
M[1,11]= sum(OffFN) #,val=True) #.shape[0]
M[1,12]= sum(RUP) #,val=True) #.shape[0]
M[1,13]= sum(RUT) #,val=True) #.shape[0]
M[1,14]= sum(RUTP) #,val=True) #.shape[0]
M[1,15]= sum(RUFP) #,val=True) #.shape[0]
M[1,16]= sum(RUTN) #,val=True) #.shape[0]
M[1,17]= sum(RUFN) #,val=True) #.shape[0]
M[1,18]= sum(RDP) #,val=True) #.shape[0]
M[1,19]= sum(RDT) #,val=True) #.shape[0]
M[1,20]= sum(RDTP) #,val=True) #.shape[0]
M[1,21]= sum(RDFP) #,val=True) #.shape[0]
M[1,22]= sum(RDTN) #,val=True) #.shape[0]
M[1,23]= sum(RDFN) #,val=True) #.shape[0]
M[1,24]= sum(OtherP) #,val=True) #.shape[0]
M[1,25]= sum(OtherT) #,val=True) #shape[0]

pd.DataFrame(M).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Count_TLs.csv')


