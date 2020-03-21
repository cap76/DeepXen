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

#Load in the raw data (labels are encoded in column 4) and one-hot encode.
Output2 = genfromtxt('allpeaks_labelled_cum_final.se.bed',delimiter='\t',dtype=None)

#Load in the list of genes
onlist = genfromtxt('OnMemFDRtGenes.txt',delimiter='\t',dtype=None)
downlist = genfromtxt('ReprogrammedDowntGenes.txt',delimiter='\t',dtype=None)
offlist = genfromtxt('OffMemFDRtGenes.txt',delimiter='\t',dtype=None)
uplist = genfromtxt('ReprogrammedUptGenes.txt',delimiter='\t',dtype=None)
complist = genfromtxt('ComplementarySettGenes.txt',delimiter='\t',dtype=None)
#foxlist = genfromtxt('FOXA1ChIP.txt',delimiter='\t',dtype=None)
#enhlist = genfromtxt('enhancers.bed',delimiter='\t',dtype=None)


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

explist = genfromtxt('edger_de_pairing_BIGtable_endoIVF_ectoDonor_endoNT_plus_DE.p.csv',delimiter='\t',dtype=None)
explist1 = genfromtxt('edger_de_pairing_BIGtable_endoIVF_ectoDonor_endoNT_plus_DE.rmnan.csv',delimiter='\t',dtype=None)

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

#¢fil=sorted(glob.glob('SomiteEcto/*.tab'))
fil=sorted(glob.glob('/Users/christopherpenfold/Desktop/Code/deepexplain/TSS_deeper/Endoderm/*.tab'))

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

from keras.models import load_model
model = load_model('Model_deeper.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model2 = load_model('Model_deeper_neg.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model3 = load_model('Model_deeper_pos.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)

model4 = load_model('Model_deeper_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model5 = load_model('Model_deeper_neg_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)
model6 = load_model('Model_deeper_pos_seed=2.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)


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
plt.savefig('PR.pdf')
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
plt.savefig('AUC.pdf')
plt.draw()
