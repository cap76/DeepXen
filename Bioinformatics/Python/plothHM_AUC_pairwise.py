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

#Normalise the expression data
Xexp = copy.deepcopy(Xexpa)
Xexp[:,0,0] = ( Xexp[:,0,0] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0])
Xexp[:,0,1] = ( Xexp[:,0,1] - np.nanmean(Xexp[trainset,0,1]) ) / np.nanstd(Xexp[trainset,0,1])

Xchra = copy.deepcopy(Xchr)

#Normalise histone modification data
for k in range(0, np.shape(fil)[0]):
      Xchr[:,0:600,k] = ( Xchr[:,0:600,k] - np.nanmean(Xchr[trainset,0:600,k]) ) / np.nanstd(Xchr[trainset,0:600,k])

plt.switch_backend('agg')

model = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model2 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_neg.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model3 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_pos.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model4 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model5 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_neg_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model6 = load_model('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Model_deeper_pos_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)

Scores = np.zeros((6,1))
Scores[0,0] = model.evaluate([Xchr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
Scores[1,0] = model3.evaluate([Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
Scores[2,0] = model2.evaluate([Xchr[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
Scores[3,0] = model4.evaluate([Xchr[valset,:,:],Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
Scores[4,0] = model6.evaluate([Xexp[valset,:,:]], Yalt[valset,:], batch_size=32)[1]
Scores[5,0] = model5.evaluate([Xchr[valset,:,:]], Yalt[valset,:], batch_size=32)[1]

pd.DataFrame(Scores, columns=['Accuracy']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/PredictionScores.csv')

#Make some predictions
classpred1 = model.predict([Xchr,Xexp], batch_size=1000)
classpred2 = model2.predict([Xchr], batch_size=1000)
classpred3 = model3.predict([Xexp], batch_size=1000)
classpred4 = model4.predict([Xchr,Xexp], batch_size=1000)
classpred5 = model5.predict([Xchr], batch_size=1000)
classpred6 = model6.predict([Xexp], batch_size=1000)

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

df14a = pd.DataFrame(classpred1,columns=['pC1','pC2','pC3','pC4','pC5'])
df14b = pd.DataFrame(classpred2,columns=['pC1','pC2','pC3','pC4','pC5'])
df14c = pd.DataFrame(classpred3,columns=['pC1','pC2','pC3','pC4','pC5'])
df14d = pd.DataFrame(classpred4,columns=['pC1','pC2','pC3','pC4','pC5'])
df14e = pd.DataFrame(classpred5,columns=['pC1','pC2','pC3','pC4','pC5'])
df14f = pd.DataFrame(classpred6,columns=['pC1','pC2','pC3','pC4','pC5'])

prediction1 = np.concatenate((df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14,df15,df16, df14a),1)
prediction2 = np.concatenate((df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14,df15,df16, df14b),1)
prediction3 = np.concatenate((df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14,df15,df16, df14c),1)
prediction4 = np.concatenate((df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14,df15,df16, df14d),1)
prediction5 = np.concatenate((df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14,df15,df16, df14e),1)
prediction6 = np.concatenate((df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14,df15,df16, df14f),1)

pd.DataFrame(prediction1, columns=['Chr','start','end','Gene','IVF','Donor','NT','FC','pVal','FC','pVal','ON','RepDown','Off','RepUp','Neg','pC1','pC2','pC3','pC4','pC5']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Pred1.csv')
pd.DataFrame(prediction2, columns=['Chr','start','end','Gene','IVF','Donor','NT','FC','pVal','FC','pVal','ON','RepDown','Off','RepUp','Neg','pC1','pC2','pC3','pC4','pC5']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Pred2Neg.csv')
pd.DataFrame(prediction3, columns=['Chr','start','end','Gene','IVF','Donor','NT','FC','pVal','FC','pVal','ON','RepDown','Off','RepUp','Neg','pC1','pC2','pC3','pC4','pC5']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Pred3Pos.csv')
pd.DataFrame(prediction4, columns=['Chr','start','end','Gene','IVF','Donor','NT','FC','pVal','FC','pVal','ON','RepDown','Off','RepUp','Neg','pC1','pC2','pC3','pC4','pC5']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Pred4.csv')
pd.DataFrame(prediction5, columns=['Chr','start','end','Gene','IVF','Donor','NT','FC','pVal','FC','pVal','ON','RepDown','Off','RepUp','Neg','pC1','pC2','pC3','pC4','pC5']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Pred5Neg.csv')
pd.DataFrame(prediction6, columns=['Chr','start','end','Gene','IVF','Donor','NT','FC','pVal','FC','pVal','ON','RepDown','Off','RepUp','Neg','pC1','pC2','pC3','pC4','pC5']).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Pred6Pos.csv')

#NOW GENERATE SOME EXAMPLE PLOTS
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
plt.switch_backend('agg')

from sklearn.metrics import roc_curve, precision_recall_curve,auc
#fpr_keras, tpr_keras, thresholds_keras = roc_curve(Yalt.ravel(), model.predict([Xchr,Xexp], batch_size=1000).ravel())
#pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt.ravel(), model.predict([Xchr,Xexp], batch_size=1000).ravel())

fpr_keras0, tpr_keras0, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_keras0, re_keras0, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())

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

np.save('EndoFull/fpr_kearst.npy', fpr_keras)
np.save('EndoFull/fpr_kears1.npy', fpr_keras1)
np.save('EndoFull/fpr_kears2.npy', fpr_keras2)
np.save('EndoFull/fpr_kears4.npy', fpr_kerast4)
np.save('EndoFull/fpr_kears5.npy', fpr_keras5)
np.save('EndoFull/fpr_kears6.npy', fpr_keras6)

np.save('EndoFull/tpr_kearst.npy', tpr_keras)
np.save('EndoFull/tpr_kears1.npy', tpr_keras1)
np.save('EndoFull/tpr_kears2.npy', tpr_keras2)
np.save('EndoFull/tpr_kears4.npy', tpr_kerast4)
np.save('EndoFull/tpr_kears5.npy', tpr_keras5)
np.save('EndoFull/tpr_kears6.npy', tpr_keras6)

np.save('EndoFull/pr_kearst.npy', pr_keras)
np.save('EndoFull/pr_kears1.npy', pr_keras1)
np.save('EndoFull/pr_kears2.npy', pr_keras2)
np.save('EndoFull/pr_kears4.npy', pr_kerast4)
np.save('EndoFull/pr_kears5.npy', pr_keras5)
np.save('EndoFull/pr_kears6.npy', pr_keras6)

np.save('EndoFull/re_kearst.npy', re_keras)
np.save('EndoFull/re_kears1.npy', re_keras1)
np.save('EndoFull/re_kears2.npy', re_keras2)
np.save('EndoFull/re_kears4.npy', re_kerast4)
np.save('EndoFull/re_kears5.npy', re_keras5)
np.save('EndoFull/re_kears6.npy', re_keras6)

from sklearn.metrics import roc_curve, precision_recall_curve,auc
#fpr_keras, tpr_keras, thresholds_keras = roc_curve(Yalt.ravel(), model.predict([Xchr,Xexp], batch_size=1000).ravel())
#pr_keras, re_keras, thresholds_keras = precision_recall_curve(Yalt.ravel(), model.predict([Xchr,Xexp], batch_size=1000).ravel())
Y0 = Yalt[valset,:]
P0 = model.predict([Xchr[valset,:,:],Xexp[valset,:,:]])
Y1 = Yalt[valset,:]
P1 = model2.predict(Xchr[valset,:,:])
Y2 = Yalt[valset,:]
P2 = model3.predict(Xexp[valset,:,:])
Y4 = Yalt[valset,:]
P4 = model4.predict([Xchr[valset,:,:],Xexp[valset,:,:]])
Y5 = Yalt[valset,:]
P5 = model5.predict(Xchr[valset,:,:])
Y6 = Yalt[valset,:]
P6 = model6.predict(Xexp[valset,:,:])

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

fpr_keras4_1, tpr_keras4_1, thresholds_kerast = roc_curve(Y4[:,[True,False,False,False,False]], P4[:,[True,False,False,False,False]])
fpr_keras4_2, tpr_keras4_2, thresholds_kerast = roc_curve(Y4[:,[False,True,False,False,False]], P4[:,[False,True,False,False,False]])
fpr_keras4_3, tpr_keras4_3, thresholds_kerast = roc_curve(Y4[:,[False,False,True,False,False]], P4[:,[False,False,True,False,False]])
fpr_keras4_4, tpr_keras4_4, thresholds_kerast = roc_curve(Y4[:,[False,False,False,True,False]], P4[:,[False,False,False,True,False]])
fpr_keras4_5, tpr_keras4_5, thresholds_kerast = roc_curve(Y4[:,[False,False,False,False,True]], P4[:,[False,False,False,False,True]])

fpr_keras5_1, tpr_keras5_1, thresholds_kerast = roc_curve(Y5[:,[True,False,False,False,False]], P5[:,[True,False,False,False,False]])
fpr_keras5_2, tpr_keras5_2, thresholds_kerast = roc_curve(Y5[:,[False,True,False,False,False]], P5[:,[False,True,False,False,False]])
fpr_keras5_3, tpr_keras5_3, thresholds_kerast = roc_curve(Y5[:,[False,False,True,False,False]], P5[:,[False,False,True,False,False]])
fpr_keras5_4, tpr_keras5_4, thresholds_kerast = roc_curve(Y5[:,[False,False,False,True,False]], P5[:,[False,False,False,True,False]])
fpr_keras5_5, tpr_keras5_5, thresholds_kerast = roc_curve(Y5[:,[False,False,False,False,True]], P5[:,[False,False,False,False,True]])

fpr_keras6_1, tpr_keras6_1, thresholds_kerast = roc_curve(Y6[:,[True,False,False,False,False]], P6[:,[True,False,False,False,False]])
fpr_keras6_2, tpr_keras6_2, thresholds_kerast = roc_curve(Y6[:,[False,True,False,False,False]], P6[:,[False,True,False,False,False]])
fpr_keras6_3, tpr_keras6_3, thresholds_kerast = roc_curve(Y6[:,[False,False,True,False,False]], P6[:,[False,False,True,False,False]])
fpr_keras6_4, tpr_keras6_4, thresholds_kerast = roc_curve(Y6[:,[False,False,False,True,False]], P6[:,[False,False,False,True,False]])
fpr_keras6_5, tpr_keras6_5, thresholds_kerast = roc_curve(Y6[:,[False,False,False,False,True]], P6[:,[False,False,False,False,True]])

np.save('EndoFull/tpr_kears0_1.npy', tpr_keras0_1)
np.save('EndoFull/tpr_kears0_2.npy', tpr_keras0_2)
np.save('EndoFull/tpr_kears0_3.npy', tpr_keras0_3)
np.save('EndoFull/tpr_kears0_4.npy', tpr_keras0_4)
np.save('EndoFull/tpr_kears0_5.npy', tpr_keras0_5)
np.save('EndoFull/tpr_kears1_1.npy', tpr_keras1_1)
np.save('EndoFull/tpr_kears1_2.npy', tpr_keras1_2)
np.save('EndoFull/tpr_kears1_3.npy', tpr_keras1_3)
np.save('EndoFull/tpr_kears1_4.npy', tpr_keras1_4)
np.save('EndoFull/tpr_kears1_5.npy', tpr_keras1_5)
np.save('EndoFull/tpr_kears2_1.npy', tpr_keras2_1)
np.save('EndoFull/tpr_kears2_2.npy', tpr_keras2_2)
np.save('EndoFull/tpr_kears2_3.npy', tpr_keras2_3)
np.save('EndoFull/tpr_kears2_4.npy', tpr_keras2_4)
np.save('EndoFull/tpr_kears2_5.npy', tpr_keras2_5)
np.save('EndoFull/tpr_kears4_1.npy', tpr_keras4_1)
np.save('EndoFull/tpr_kears4_2.npy', tpr_keras4_2)
np.save('EndoFull/tpr_kears4_3.npy', tpr_keras4_3)
np.save('EndoFull/tpr_kears4_4.npy', tpr_keras4_4)
np.save('EndoFull/tpr_kears4_5.npy', tpr_keras4_5)
np.save('EndoFull/tpr_kears5_1.npy', tpr_keras5_1)
np.save('EndoFull/tpr_kears5_2.npy', tpr_keras5_2)
np.save('EndoFull/tpr_kears5_3.npy', tpr_keras5_3)
np.save('EndoFull/tpr_kears5_4.npy', tpr_keras5_4)
np.save('EndoFull/tpr_kears5_5.npy', tpr_keras5_5)
np.save('EndoFull/tpr_kears6_1.npy', tpr_keras6_1)
np.save('EndoFull/tpr_kears6_2.npy', tpr_keras6_2)
np.save('EndoFull/tpr_kears6_3.npy', tpr_keras6_3)
np.save('EndoFull/tpr_kears6_4.npy', tpr_keras6_4)
np.save('EndoFull/tpr_kears6_5.npy', tpr_keras6_5)


np.save('EndoFull/fpr_kears0_1.npy', fpr_keras0_1)
np.save('EndoFull/fpr_kears0_2.npy', fpr_keras0_2)
np.save('EndoFull/fpr_kears0_3.npy', fpr_keras0_3)
np.save('EndoFull/fpr_kears0_4.npy', fpr_keras0_4)
np.save('EndoFull/fpr_kears0_5.npy', fpr_keras0_5)
np.save('EndoFull/fpr_kears1_1.npy', fpr_keras1_1)
np.save('EndoFull/fpr_kears1_2.npy', fpr_keras1_2)
np.save('EndoFull/fpr_kears1_3.npy', fpr_keras1_3)
np.save('EndoFull/fpr_kears1_4.npy', fpr_keras1_4)
np.save('EndoFull/fpr_kears1_5.npy', fpr_keras1_5)
np.save('EndoFull/fpr_kears2_1.npy', fpr_keras2_1)
np.save('EndoFull/fpr_kears2_2.npy', fpr_keras2_2)
np.save('EndoFull/fpr_kears2_3.npy', fpr_keras2_3)
np.save('EndoFull/fpr_kears2_4.npy', fpr_keras2_4)
np.save('EndoFull/fpr_kears2_5.npy', fpr_keras2_5)
np.save('EndoFull/fpr_kears4_1.npy', fpr_keras4_1)
np.save('EndoFull/fpr_kears4_2.npy', fpr_keras4_2)
np.save('EndoFull/fpr_kears4_3.npy', fpr_keras4_3)
np.save('EndoFull/fpr_kears4_4.npy', fpr_keras4_4)
np.save('EndoFull/fpr_kears4_5.npy', fpr_keras4_5)
np.save('EndoFull/fpr_kears5_1.npy', fpr_keras5_1)
np.save('EndoFull/fpr_kears5_2.npy', fpr_keras5_2)
np.save('EndoFull/fpr_kears5_3.npy', fpr_keras5_3)
np.save('EndoFull/fpr_kears5_4.npy', fpr_keras5_4)
np.save('EndoFull/fpr_kears5_5.npy', fpr_keras5_5)
np.save('EndoFull/fpr_kears6_1.npy', fpr_keras6_1)
np.save('EndoFull/fpr_kears6_2.npy', fpr_keras6_2)
np.save('EndoFull/fpr_kears6_3.npy', fpr_keras6_3)
np.save('EndoFull/fpr_kears6_4.npy', fpr_keras6_4)
np.save('EndoFull/fpr_kears6_5.npy', fpr_keras6_5)


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

pr_keras4_1, re_keras4_1, thresholds_kerast = precision_recall_curve(Y4[:,[True,False,False,False,False]], P4[:,[True,False,False,False,False]])
pr_keras4_2, re_keras4_2, thresholds_kerast = precision_recall_curve(Y4[:,[False,True,False,False,False]], P4[:,[False,True,False,False,False]])
pr_keras4_3, re_keras4_3, thresholds_kerast = precision_recall_curve(Y4[:,[False,False,True,False,False]], P4[:,[False,False,True,False,False]])
pr_keras4_4, re_keras4_4, thresholds_kerast = precision_recall_curve(Y4[:,[False,False,False,True,False]], P4[:,[False,False,False,True,False]])
pr_keras4_5, re_keras4_5, thresholds_kerast = precision_recall_curve(Y4[:,[False,False,False,False,True]], P4[:,[False,False,False,False,True]])

pr_keras5_1, re_keras5_1, thresholds_kerast = precision_recall_curve(Y5[:,[True,False,False,False,False]], P5[:,[True,False,False,False,False]])
pr_keras5_2, re_keras5_2, thresholds_kerast = precision_recall_curve(Y5[:,[False,True,False,False,False]], P5[:,[False,True,False,False,False]])
pr_keras5_3, re_keras5_3, thresholds_kerast = precision_recall_curve(Y5[:,[False,False,True,False,False]], P5[:,[False,False,True,False,False]])
pr_keras5_4, re_keras5_4, thresholds_kerast = precision_recall_curve(Y5[:,[False,False,False,True,False]], P5[:,[False,False,False,True,False]])
pr_keras5_5, re_keras5_5, thresholds_kerast = precision_recall_curve(Y5[:,[False,False,False,False,True]], P5[:,[False,False,False,False,True]])

pr_keras6_1, re_keras6_1, thresholds_kerast = precision_recall_curve(Y6[:,[True,False,False,False,False]], P6[:,[True,False,False,False,False]])
pr_keras6_2, re_keras6_2, thresholds_kerast = precision_recall_curve(Y6[:,[False,True,False,False,False]], P6[:,[False,True,False,False,False]])
pr_keras6_3, re_keras6_3, thresholds_kerast = precision_recall_curve(Y6[:,[False,False,True,False,False]], P6[:,[False,False,True,False,False]])
pr_keras6_4, re_keras6_4, thresholds_kerast = precision_recall_curve(Y6[:,[False,False,False,True,False]], P6[:,[False,False,False,True,False]])
pr_keras6_5, re_keras6_5, thresholds_kerast = precision_recall_curve(Y6[:,[False,False,False,False,True]], P6[:,[False,False,False,False,True]])

np.save('EndoFull/pr_kears0_1.npy', pr_keras0_1)
np.save('EndoFull/pr_kears0_2.npy', pr_keras0_2)
np.save('EndoFull/pr_kears0_3.npy', pr_keras0_3)
np.save('EndoFull/pr_kears0_4.npy', pr_keras0_4)
np.save('EndoFull/pr_kears0_5.npy', pr_keras0_5)
np.save('EndoFull/pr_kears1_1.npy', pr_keras1_1)
np.save('EndoFull/pr_kears1_2.npy', pr_keras1_2)
np.save('EndoFull/pr_kears1_3.npy', pr_keras1_3)
np.save('EndoFull/pr_kears1_4.npy', pr_keras1_4)
np.save('EndoFull/pr_kears1_5.npy', pr_keras1_5)
np.save('EndoFull/pr_kears2_1.npy', pr_keras2_1)
np.save('EndoFull/pr_kears2_2.npy', pr_keras2_2)
np.save('EndoFull/pr_kears2_3.npy', pr_keras2_3)
np.save('EndoFull/pr_kears2_4.npy', pr_keras2_4)
np.save('EndoFull/pr_kears2_5.npy', pr_keras2_5)
np.save('EndoFull/pr_kears4_1.npy', pr_keras4_1)
np.save('EndoFull/pr_kears4_2.npy', pr_keras4_2)
np.save('EndoFull/pr_kears4_3.npy', pr_keras4_3)
np.save('EndoFull/pr_kears4_4.npy', pr_keras4_4)
np.save('EndoFull/pr_kears4_5.npy', pr_keras4_5)
np.save('EndoFull/pr_kears5_1.npy', pr_keras5_1)
np.save('EndoFull/pr_kears5_2.npy', pr_keras5_2)
np.save('EndoFull/pr_kears5_3.npy', pr_keras5_3)
np.save('EndoFull/pr_kears5_4.npy', pr_keras5_4)
np.save('EndoFull/pr_kears5_5.npy', pr_keras5_5)
np.save('EndoFull/pr_kears6_1.npy', pr_keras6_1)
np.save('EndoFull/pr_kears6_2.npy', pr_keras6_2)
np.save('EndoFull/pr_kears6_3.npy', pr_keras6_3)
np.save('EndoFull/pr_kears6_4.npy', pr_keras6_4)
np.save('EndoFull/pr_kears6_5.npy', pr_keras6_5)

np.save('EndoFull/re_kears0_1.npy', re_keras0_1)
np.save('EndoFull/re_kears0_2.npy', re_keras0_2)
np.save('EndoFull/re_kears0_3.npy', re_keras0_3)
np.save('EndoFull/re_kears0_4.npy', re_keras0_4)
np.save('EndoFull/re_kears0_5.npy', re_keras0_5)
np.save('EndoFull/re_kears1_1.npy', re_keras1_1)
np.save('EndoFull/re_kears1_2.npy', re_keras1_2)
np.save('EndoFull/re_kears1_3.npy', re_keras1_3)
np.save('EndoFull/re_kears1_4.npy', re_keras1_4)
np.save('EndoFull/re_kears1_5.npy', re_keras1_5)
np.save('EndoFull/re_kears2_1.npy', re_keras2_1)
np.save('EndoFull/re_kears2_2.npy', re_keras2_2)
np.save('EndoFull/re_kears2_3.npy', re_keras2_3)
np.save('EndoFull/re_kears2_4.npy', re_keras2_4)
np.save('EndoFull/re_kears2_5.npy', re_keras2_5)
np.save('EndoFull/re_kears4_1.npy', re_keras4_1)
np.save('EndoFull/re_kears4_2.npy', re_keras4_2)
np.save('EndoFull/re_kears4_3.npy', re_keras4_3)
np.save('EndoFull/re_kears4_4.npy', re_keras4_4)
np.save('EndoFull/re_kears4_5.npy', re_keras4_5)
np.save('EndoFull/re_kears5_1.npy', re_keras5_1)
np.save('EndoFull/re_kears5_2.npy', re_keras5_2)
np.save('EndoFull/re_kears5_3.npy', re_keras5_3)
np.save('EndoFull/re_kears5_4.npy', re_keras5_4)
np.save('EndoFull/re_kears5_5.npy', re_keras5_5)
np.save('EndoFull/re_kears6_1.npy', re_keras6_1)
np.save('EndoFull/re_kears6_2.npy', re_keras6_2)
np.save('EndoFull/re_kears6_3.npy', re_keras6_3)
np.save('EndoFull/re_kears6_4.npy', re_keras6_4)
np.save('EndoFull/re_kears6_5.npy', re_keras6_5)

#Y0_1
#Y0_2
#Y0_3
#Y0_4
#Y0_5






fpr_keras0, tpr_keras0, thresholds_kerast = roc_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())
pr_keras0, re_keras0, thresholds_keras = precision_recall_curve(Yalt[valset,:].ravel(), model.predict([Xchr[valset,:,:],Xexp[valset,:,:]], batch_size=1000).ravel())

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
plt.plot(re_keras1, pr_keras1, label='Chromatin') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(re_keras2, pr_keras2, label='Expression') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(re_kerast4, pr_kerast4, label='Combined (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(re_keras5, pr_keras5, label='Chromatin (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
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
plt.plot(fpr_keras1, tpr_keras1, label='Chromatin') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(fpr_keras2, tpr_keras2, label='Expression') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(fpr_kerast4, tpr_kerast4, label='Combined (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
plt.plot(fpr_keras5, tpr_keras5, label='Chromatin (2)') #, label='ROC curve (area = %0.2f)' % roc_auc[i])
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
classpred = classpred3
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

M[2,0]= sum(OnP) #,val=True)
M[2,1]= sum(OnT) #,val=True) #.shape[0]
M[2,2]= sum(OnTP) #,val=True) #.shape[0]
M[2,3]= sum(OnFP) #,val=True) #.shape[0]
M[2,4]= sum(OnTN) #,val=True) #.shape[0]
M[2,5]= sum(OnFN) #,val=True) #.shape[0]
M[2,6]= sum(OffP) #,val=True) #.shape[0]
M[2,7]= sum(OffT) #,val=True) #.shape[0]
M[2,8]= sum(OffTP) #,val=True) #.shape[0]
M[2,9]= sum(OffFP) #,val=True) #.shape[0]
M[2,10]= sum(OnTN) #,val=True) #.shape[0]
M[2,11]= sum(OffFN) #,val=True) #.shape[0]
M[2,12]= sum(RUP) #,val=True) #.shape[0]
M[2,13]= sum(RUT) #,val=True) #.shape[0]
M[2,14]= sum(RUTP) #,val=True) #.shape[0]
M[2,15]= sum(RUFP) #,val=True) #.shape[0]
M[2,16]= sum(RUTN) #,val=True) #.shape[0]
M[2,17]= sum(RUFN) #,val=True) #.shape[0]
M[2,18]= sum(RDP) #,val=True) #.shape[0]
M[2,19]= sum(RDT) #,val=True) #.shape[0]
M[2,20]= sum(RDTP) #,val=True) #.shape[0]
M[2,21]= sum(RDFP) #,val=True) #.shape[0]
M[2,22]= sum(RDTN) #,val=True) #.shape[0]
M[2,23]= sum(RDFN) #,val=True) #.shape[0]
M[2,24]= sum(OtherP) #,val=True) #.shape[0]
M[2,25]= sum(OtherT) #,val=True) #shape[0]

classpred = classpred4
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

M[3,0]= sum(OnP) #,val=True)
M[3,1]= sum(OnT) #,val=True) #.shape[0]
M[3,2]= sum(OnTP) #,val=True) #.shape[0]
M[3,3]= sum(OnFP) #,val=True) #.shape[0]
M[3,4]= sum(OnTN) #,val=True) #.shape[0]
M[3,5]= sum(OnFN) #,val=True) #.shape[0]
M[3,6]= sum(OffP) #,val=True) #.shape[0]
M[3,7]= sum(OffT) #,val=True) #.shape[0]
M[3,8]= sum(OffTP) #,val=True) #.shape[0]
M[3,9]= sum(OffFP) #,val=True) #.shape[0]
M[3,10]= sum(OnTN) #,val=True) #.shape[0]
M[3,11]= sum(OffFN) #,val=True) #.shape[0]
M[3,12]= sum(RUP) #,val=True) #.shape[0]
M[3,13]= sum(RUT) #,val=True) #.shape[0]
M[3,14]= sum(RUTP) #,val=True) #.shape[0]
M[3,15]= sum(RUFP) #,val=True) #.shape[0]
M[3,16]= sum(RUTN) #,val=True) #.shape[0]
M[3,17]= sum(RUFN) #,val=True) #.shape[0]
M[3,18]= sum(RDP) #,val=True) #.shape[0]
M[3,19]= sum(RDT) #,val=True) #.shape[0]
M[3,20]= sum(RDTP) #,val=True) #.shape[0]
M[3,21]= sum(RDFP) #,val=True) #.shape[0]
M[3,22]= sum(RDTN) #,val=True) #.shape[0]
M[3,23]= sum(RDFN) #,val=True) #.shape[0]
M[3,24]= sum(OtherP) #,val=True) #.shape[0]
M[3,25]= sum(OtherT) #,val=True) #shape[0]



classpred = classpred5
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


M[4,0]= sum(OnP) #,val=True)
M[4,1]= sum(OnT) #,val=True) #.shape[0]
M[4,2]= sum(OnTP) #,val=True) #.shape[0]
M[4,3]= sum(OnFP) #,val=True) #.shape[0]
M[4,4]= sum(OnTN) #,val=True) #.shape[0]
M[4,5]= sum(OnFN) #,val=True) #.shape[0]
M[4,6]= sum(OffP) #,val=True) #.shape[0]
M[4,7]= sum(OffT) #,val=True) #.shape[0]
M[4,8]= sum(OffTP) #,val=True) #.shape[0]
M[4,9]= sum(OffFP) #,val=True) #.shape[0]
M[4,10]= sum(OnTN) #,val=True) #.shape[0]
M[4,11]= sum(OffFN) #,val=True) #.shape[0]
M[4,12]= sum(RUP) #,val=True) #.shape[0]
M[4,13]= sum(RUT) #,val=True) #.shape[0]
M[4,14]= sum(RUTP) #,val=True) #.shape[0]
M[4,15]= sum(RUFP) #,val=True) #.shape[0]
M[4,16]= sum(RUTN) #,val=True) #.shape[0]
M[4,17]= sum(RUFN) #,val=True) #.shape[0]
M[4,18]= sum(RDP) #,val=True) #.shape[0]
M[4,19]= sum(RDT) #,val=True) #.shape[0]
M[4,20]= sum(RDTP) #,val=True) #.shape[0]
M[4,21]= sum(RDFP) #,val=True) #.shape[0]
M[4,22]= sum(RDTN) #,val=True) #.shape[0]
M[4,23]= sum(RDFN) #,val=True) #.shape[0]
M[4,24]= sum(OtherP) #,val=True) #.shape[0]
M[4,25]= sum(OtherT) #,val=True) #shape[0]

classpred = classpred6
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


M[5,0]= sum(OnP) #,val=True)
M[5,1]= sum(OnT) #,val=True) #.shape[0]
M[5,2]= sum(OnTP) #,val=True) #.shape[0]
M[5,3]= sum(OnFP) #,val=True) #.shape[0]
M[5,4]= sum(OnTN) #,val=True) #.shape[0]
M[5,5]= sum(OnFN) #,val=True) #.shape[0]
M[5,6]= sum(OffP) #,val=True) #.shape[0]
M[5,7]= sum(OffT) #,val=True) #.shape[0]
M[5,8]= sum(OffTP) #,val=True) #.shape[0]
M[5,9]= sum(OffFP) #,val=True) #.shape[0]
M[5,10]= sum(OnTN) #,val=True) #.shape[0]
M[5,11]= sum(OffFN) #,val=True) #.shape[0]
M[5,12]= sum(RUP) #,val=True) #.shape[0]
M[5,13]= sum(RUT) #,val=True) #.shape[0]
M[5,14]= sum(RUTP) #,val=True) #.shape[0]
M[5,15]= sum(RUFP) #,val=True) #.shape[0]
M[5,16]= sum(RUTN) #,val=True) #.shape[0]
M[5,17]= sum(RUFN) #,val=True) #.shape[0]
M[5,18]= sum(RDP) #,val=True) #.shape[0]
M[5,19]= sum(RDT) #,val=True) #.shape[0]
M[5,20]= sum(RDTP) #,val=True) #.shape[0]
M[5,21]= sum(RDFP) #,val=True) #.shape[0]
M[5,22]= sum(RDTN) #,val=True) #.shape[0]
M[5,23]= sum(RDFN) #,val=True) #.shape[0]
M[5,24]= sum(OtherP) #,val=True) #.shape[0]
M[5,25]= sum(OtherT) #,val=True) #shape[0]

pd.DataFrame(M).to_csv('/mnt/scratch/gurdon/cap76/DeepXen/ResultsEndoAll/Counts.csv')

