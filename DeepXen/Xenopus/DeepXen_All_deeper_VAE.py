# Multiple Inputs
import keras
from keras import layers

#from keras.utils import plot_model
#from keras.models import Model
#from keras.layers import Input
#from keras.layers import Dense
#from keras.layers import Flatten
#from keras.layers import Activation
#from keras.layers import Dropout
#from keras import optimizers
#from keras import regularizers
#from keras.layers.convolutional import Conv2D
#from keras.layers.convolutional import Conv1D
#from keras.layers.pooling import MaxPooling2D
#from keras.layers.pooling import MaxPooling1D
#from keras.layers.merge import concatenate

#from keras.utils import plot_model
#from keras.models import Model
#from keras.layers import Input
#from keras.layers import Dense
#from keras.layers import Flatten
#from keras.layers import Activation
#from keras.layers import Dropout
#from keras import optimizers
#from keras import regularizers
#from keras.layers.convolutional import Conv2D
#from keras.layers.convolutional import Conv1D
#from keras.layers.pooling import MaxPooling2D
#from keras.layers.pooling import MaxPooling1D
#from keras.layers.merge import concatenate
from keras.utils import np_utils
from sklearn.preprocessing import LabelEncoder
from numpy import genfromtxt
import numpy as np
import pandas as pd
import glob
import os
import sys
import copy
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
Output2 = genfromtxt('BW_Endo_All2/allpeaks_labelled_cum_final.se.bed',delimiter='\t',dtype=None)
#Load in the list of genes
onlist = genfromtxt('ChIP_Endo_All/OnMemFDRtGenes.txt',delimiter='\t',dtype=None)
downlist = genfromtxt('ChIP_Endo_All/ReprogrammedDowntGenes.txt',delimiter='\t',dtype=None)
offlist = genfromtxt('ChIP_Endo_All/OffMemFDRtGenes.txt',delimiter='\t',dtype=None)
uplist = genfromtxt('ChIP_Endo_All/ReprogrammedUptGenes.txt',delimiter='\t',dtype=None)
complist = genfromtxt('ChIP_Endo_All/ComplementarySettGenes.txt',delimiter='\t',dtype=None)
#Generate the Y matrix
Yalt = np.zeros((np.shape(Output2)[0],5))
Xexpa = np.zeros((np.shape(Output2)[0],1,3))

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
   Xexpa[i,0,2] = np.log2(explist['f3'][arr_ind].astype('f')+1)

np.nan_to_num(Xexpa,copy=False)
fil1=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/Endo*bwe.tab'))
fil2=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/Ecto*bwe.tab'))
fil3=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/GC*bwe.tab'))
fil4=sorted(glob.glob('/mnt/scratch/gurdon/cap76/DeepXen/BW_Endo_All2/Meth*bwe.tab'))
fil_1 = fil1 + fil3 # + fil4
fil_2 = fil2 + fil3 # + fil4

#meh = 0
#for k in range(0, np.shape(fil_1)[0]):
#      meh = meh +1
#Load in all the expression data
encoder = LabelEncoder()
meh = encoder.fit(Output2['f0'])
Xchr = np.zeros((np.shape(Output2)[0],600,np.shape(fil_1)[0]+2))
for k in range(0, np.shape(fil_1)[0]):
        Output = genfromtxt(fil_1[k],delimiter='\t',dtype=None, skip_header=3)
        Xchr[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output
kadd = k

#Load in all the expression data
encoder = LabelEncoder()
meh = encoder.fit(Output2['f0'])
XchrO = np.zeros((np.shape(Output2)[0],600,np.shape(fil_2)[0]+2))
for k in range(0, np.shape(fil_2)[0]):
        Output = genfromtxt(fil_2[k],delimiter='\t',dtype=None, skip_header=3)
        XchrO[0:np.shape(Output)[0],0:np.shape(Output)[1],k] = Output

#Split into train/val/test sets
trainset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-18]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-17]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-16]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-15]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-14]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-13]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-12]))
valset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-11]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-10]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-9]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-8]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-7]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-6]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-5]))
testset=((Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-3]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-2]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-1]) | (Output2['f0']==meh.classes_[np.shape(meh.classes_)[0]-4]))

#Normalise the expression data
Xexp = copy.deepcopy(Xexpa)
Xexp[:,0,1] = ( Xexp[:,0,1] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0]) #To
Xexp[:,0,2] = ( Xexp[:,0,2] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0]) #NT
Xexp[:,0,0] = ( Xexp[:,0,0] - np.nanmean(Xexp[trainset,0,0]) ) / np.nanstd(Xexp[trainset,0,0]) #From

#Xchr[1,0:3,k+1]
Xchr[0:np.shape(Output)[0],:,kadd+1] = np.transpose(np.tile(Xexp[:,0,0],(600,1)) )
Xchr[0:np.shape(Output)[0],:,kadd+2] =  np.transpose(np.tile(Xexp[:,0,1],(600,1)) )
#XchrO[0:np.shape(Output)[0],:,k+1] = Xexp[:,0,0]
#XchrO[0:np.shape(Output)[0],:,k+2] = Xexp[:,0,2]

#Xchra = copy.deepcopy(Xchr)

#Normalise histone modification data
for k in range(0, np.shape(fil_2)[0]):
      XchrO[:,0:600,k] = ( XchrO[:,0:600,k] - np.nanmean(Xchr[trainset,0:600,k]) ) / np.nanstd(Xchr[trainset,0:600,k])

for k in range(0, np.shape(fil_1)[0]):
      Xchr[:,0:600,k] = ( Xchr[:,0:600,k] - np.nanmean(Xchr[trainset,0:600,k]) ) / np.nanstd(Xchr[trainset,0:600,k])

#for k in range(0, np.shape(fil_2)[0]):
#     XchrO[:,0:600,k] = ( XchrO[:,0:600,k] - np.nanmean(XchrO[trainset,0:600,k]) ) / np.nanstd(XchrO[trainset,0:600,k])
#nclasses = np.sum(Yalt[trainset,:],0)
#cwe = np.ndarray.max(nclasses) / nclasses 
#class_weight={0: cwe[0], 1: cwe[1], 2: cwe[2], 3: cwe[3], 4: cwe[4]} #, 5: cwe[5]} #, 6: cwe[6]
#import keras
#from keras import layers

original_dim = np.shape(Xchr)[1]*np.shape(Xchr)[2]
original_dimO = np.shape(XchrO)[1]*np.shape(XchrO)[2]

latent_dim = 2
intermediate_dim = 50

inputs = keras.Input(shape=(original_dim,))
h = layers.Dense(intermediate_dim, activation='relu')(inputs)
z_mean = layers.Dense(latent_dim)(h)
z_log_sigma = layers.Dense(latent_dim)(h)

#Now run a conv autoencoder to do the thing.
#input_img = keras.Input(shape=(600, 13))
#x = layers.Conv1D(64, 10, activation='relu', padding='same')(input_img)
#x = layers.MaxPooling1D(pool_size=2, padding='same')(x)
#x = layers.Conv1D(64, 10, activation='relu', padding='same')(x)
#x = layers.MaxPooling1D(pool_size=2, padding='same')(x)
#x = layers.Conv1D(64, 10, activation='relu', padding='same')(x)
#encoded = layers.MaxPooling1D(pool_size=3, padding='same')(x)

from keras import backend as K

def sampling(args):
    z_mean, z_log_sigma = args
    epsilon = K.random_normal(shape=(K.shape(z_mean)[0], latent_dim),
                              mean=0., stddev=0.1)
    return z_mean + K.exp(z_log_sigma) * epsilon

z = layers.Lambda(sampling)([z_mean, z_log_sigma])

# Create encoder
encoder = keras.Model(inputs, [z_mean, z_log_sigma, z], name='encoder')

# Create decoder
latent_inputs = keras.Input(shape=(latent_dim,), name='z_sampling')
x = layers.Dense(intermediate_dim, activation='relu')(latent_inputs)
outputs = layers.Dense(original_dimO, activation='sigmoid')(x)
decoder = keras.Model(latent_inputs, outputs, name='decoder')

# instantiate VAE model
outputs = decoder(encoder(inputs)[2])
vae = keras.Model(inputs, outputs, name='vae_mlp')

reconstruction_loss = keras.losses.binary_crossentropy(inputs, outputs)
reconstruction_loss *= original_dim
kl_loss = 1 + z_log_sigma - K.square(z_mean) - K.exp(z_log_sigma)
kl_loss = K.sum(kl_loss, axis=-1)
kl_loss *= -0.5
vae_loss = K.mean(reconstruction_loss + kl_loss)
vae.add_loss(vae_loss)
vae.compile(optimizer='adam')

#(x_train, y_train), (x_test, y_test) = mnist.load_data()

#img_rows = 600
#img_cols = 14

# Compute VAE loss
#def my_vae_loss(y_true, y_pred):
#    xent_loss = img_rows * img_cols * metrics.binary_crossentropy(K.flatten(y_true), K.flatten(y_pred))
#    kl_loss = - 0.5 * K.sum(1 + z_log_var - K.square(z_mean) - K.exp(z_log_var), axis=-1)
#    vae_loss = K.mean(xent_loss + kl_loss)
#    return vae_loss

#vae.compile(optimizer='rmsprop', loss=my_vae_loss)


#x_train = x_train.astype('float32') / 255.
#x_test = x_test.astype('float32') / 255.

x_train = Xchr[trainset,:,:]
y_train = XchrO[trainset,:,:]

x_test = Xchr[testset,:,:]
y_test = XchrO[testset,:,:]



x_train = x_train.reshape((len(x_train), np.prod(x_train.shape[1:])))
y_train = y_train.reshape((len(y_train), np.prod(y_train.shape[1:])))

x_test = x_test.reshape((len(x_test), np.prod(x_test.shape[1:])))
y_test = y_test.reshape((len(y_test), np.prod(y_test.shape[1:])))


vae.fit( x_train, y_train,
        epochs=2,
        batch_size=32),
        validation_data=(x_test, y_test))


#x_train = x_train.reshape((len(x_train), np.prod(x_train.shape[1:])))
#x_test = x_test.reshape((len(x_test), np.prod(x_test.shape[1:])))

#vae.fit(x_train, x_train,
#        epochs=100,
#        batch_size=32,
#        validation_data=(x_test, x_test))



# Create encoder
#encoder = keras.Model(inputs, [z_mean, z_log_sigma, z], name='encoder')

# Create decoder
#latent_inputs = keras.Input(shape=(latent_dim,), name='z_sampling')
#x = layers.Dense(intermediate_dim, activation='relu')(latent_inputs)
#outputs = layers.Dense(original_dim, activation='sigmoid')(x)
#decoder = keras.Model(latent_inputs, outputs, name='decoder')

# instantiate VAE model
#outputs = decoder(encoder(inputs)[2])
#vae = keras.Model(inputs, outputs, name='vae_mlp')

#x = layers.Conv1D(64, 10, activation='relu', padding='same')(encoded)
#x = layers.UpSampling1D(size=3)(x)
#x = layers.Conv1D(64, 10, activation='relu', padding='same')(x)
#x = layers.UpSampling1D(size=2)(x)
#x = layers.Conv1D(64, 10, activation='relu', padding='same')(x)
#x = layers.UpSampling1D(size=2)(x)
decoded = layers.Conv1D(10, 10, activation='relu', padding='same')(x)
#decoded = layers.Conv1D(12, 10, activation='sigmoid', padding='same')(x)

autoencoder = keras.Model(input_img, decoded)
#autoencoder.compile(optimizer='adam', loss='binary_crossentropy')
autoencoder.compile(optimizer='adam', loss='mean_squared_error')
autoencoder.summary()

#checkpoint_filepath = 'CAE/'
#model_checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(
#    filepath=checkpoint_filepath,
#    save_weights_only=False,
#    monitor='val_mean_squared_error',
#    mode='min',
#    save_best_only=True)

from keras.callbacks import ModelCheckpoint

checkpoint = ModelCheckpoint('CAE/EndoCAE.h5', monitor='val_loss', verbose=1, save_best_only=True, mode='min')
callbacks_list = [checkpoint]


autoencoder.fit(Xchr[valset,:,:], XchrO[valset,:],
                epochs=500,
                batch_size=500,
                validation_data=(Xchr[testset,:,:], XchrO[testset,:]),
                callbacks=callbacks_list)


Xpredval = autoencoder.predict(Xchr[valset,:,:], batch_size=100)

np.save('CAE/Xpredval.npy', Xpredval)
np.save('CAE/Xval.npy', Xchr[valset,:,:])



#Now redo prediction for the NT
Xchr[0:np.shape(Output)[0],:,kadd+2] =  np.transpose(np.tile(Xexp[:,0,2],(600,1)) )
Xpredval = autoencoder.predict(Xchr, batch_size=100)
np.save('CAE/XpredvalNT.npy', Xpredval)



#autoencoder.compile(loss={'autoencoder': 'binary_crossentropy',
#                    'classification': 'binary_crossentropy'},
#              loss_weights={'autoencoder': 1.0,
#                            'classification': 1.0},
#              optimizer='adam',
#              metrics={'autoencoder': ['binary_crossentropy', 'mse'], 'classification': 'binary_accuracy'})

#checkpoint = ModelCheckpoint('./PredictEpiLC.h5', monitor='val_loss', verbose=1, save_best_only=True, mode='min')
#callbacks_list = [checkpoint]

#Just a trial run to predict classifiers. Thus will be converted to a regression task to predict Brd4 iBET levels
#autoencoder.fit(Xchr_noisy[trainset,:,:], 
#          {'autoencoder': Xchr1[trainset,:,:], 'classification': Yalt_1[trainset]},
#          batch_size=500,
#          epochs=500,
#          validation_data= (Xchr_noisy[testset,:,:], {'autoencoder': Xchr1[testset,:,:], 'classification': Yalt_1[testset]}),
#          callbacks=callbacks_list,verbose=1)


#Check we're using the best model
#from keras.models import load_model
#autoencoder = load_model('PredictEpiLC.h5') #predictions = model.predict([Xchr,Xexp], batch_size=1000)

#SmoothX,Classpred = autoencoder.predict(Xchr_noisy, batch_size=100)
#SmoothXp,Classpredp = autoencoder.predict(Xchr_noisy2, batch_size=100)
#SmoothXp2,Classpredp2 = autoencoder.predict(Xchr3, batch_size=100)

#np.savez('SmoothXp',SmoothXp)
#np.savez('SmoothX',SmoothX)
#np.savez('SmoothX',SmoothXp2)


#np.savez('Classpred',Classpred)
#np.savez('Classpredp',Classpredp)
#np.savez('Classpredp',Classpredp2)

