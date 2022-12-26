import tensorflow as tf
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Conv2D, Conv2DTranspose, Dense, Flatten, Dropout, BatchNormalization, Reshape, LeakyReLU
import matplotlib.pyplot as plt
import os
import pandas as pd
from sklearn.svm import SVC
import numpy as np
from tensorflow.keras.models import load_model


epoch = 1


compl_path = './compl'
compl_file_list = os.listdir(compl_path)
compl_file_list_py = [
    file for file in compl_file_list if file.endswith('.xlsx')]

heat_path = './heat'
heat_file_list = os.listdir(heat_path)
heat_file_list_py = [file for file in heat_file_list if file.endswith('.xlsx')]

valid_path = './valid'
valid_file_list = os.listdir(valid_path)
valid_file_list_py = [file for file in valid_file_list if file.endswith('.xlsx')]


X_compl = []
for i in compl_file_list_py:
    filename = compl_path+'/'+i
    df = pd.read_excel(filename, engine='openpyxl',
                       header=None, index_col=None)
    sample = df.to_numpy()
    fliped_sample = np.flip(df.to_numpy(), axis=1)
    for j in range(4):
        sample = np.rot90(sample,j)
        fliped_sample = np.rot90(fliped_sample,j)
        X_compl.append(sample)
        X_compl.append(fliped_sample)

X_heat = []
for i in heat_file_list_py:
    filename = heat_path+'/'+i
    df = pd.read_excel(filename, engine='openpyxl',
                       header=None, index_col=None)
    sample = df.to_numpy()
    fliped_sample = np.flip(df.to_numpy(), axis=1)
    for j in range(4):
        sample = np.rot90(sample,j)
        fliped_sample = np.rot90(fliped_sample,j)
        X_heat.append(sample)
        X_heat.append(fliped_sample)

X_valid = []
for i in valid_file_list_py:
    filename = valid_path+'/'+i
    df = pd.read_excel(filename, engine='openpyxl',
                       header=None, index_col=None)
    sample = df.to_numpy()
    X_valid.append(sample)

validation_set = np.stack(X_valid, axis=0)
validation_set = validation_set.reshape(-1, 100, 100, 1)
    


X_list = X_compl + X_heat
x_train = np.stack(X_list, axis=0)
x_train = x_train.reshape(-1, 100, 100, 1)

y_train = []
Y_compl = [0 for i in range(len(X_compl))]
Y_heat = [1 for i in range(len(X_heat))]
y_train = Y_compl + Y_heat
y_train = np.array(y_train)




########################### encoder ###########################
Input = tf.keras.Input
encoder_input = Input(shape=(100, 100, 1))

# 28 X 28
x = Conv2D(32, 3, padding='same')(encoder_input)
x = BatchNormalization()(x)
x = LeakyReLU()(x)

# 28 X 28 -> 14 X 14
x = Conv2D(64, 3, strides=2, padding='same')(x)
x = BatchNormalization()(x)
x = LeakyReLU()(x)

# 14 X 14 -> 7 X 7
x = Conv2D(64, 3, strides=2, padding='same')(x)
x = BatchNormalization()(x)
x = LeakyReLU()(x)

# 7 X 7
x = Conv2D(64, 3, padding='same')(x)
x = BatchNormalization()(x)
x = LeakyReLU()(x)

x = Flatten()(x)

# 2D 좌표로 표기하기 위하여 2를 출력값으로 지정합니다.
encoder_output = Dense(2)(x)

encoder = Model(encoder_input, encoder_output)
encoder.summary()
###############################################################


########################### dncoder ###########################
# Input으로는 2D 좌표가 들어갑니다.
decoder_input = Input(shape=(2, ))

# 2D 좌표를 7*7*64 개의 neuron 출력 값을 가지도록 변경합니다.
x = Dense(25*25*64)(decoder_input)
x = Reshape((25, 25, 64))(x)

# 7 X 7 -> 7 X 7
x = Conv2DTranspose(64, 3, strides=1, padding='same')(x)
x = BatchNormalization()(x)
x = LeakyReLU()(x)

# 7 X 7 -> 14 X 14
x = Conv2DTranspose(64, 3, strides=2, padding='same')(x)
x = BatchNormalization()(x)
x = LeakyReLU()(x)

# 14 X 14 -> 28 X 28
x = Conv2DTranspose(64, 3, strides=2, padding='same')(x)
x = BatchNormalization()(x)
x = LeakyReLU()(x)

# 28 X 28 -> 28 X 28
x = Conv2DTranspose(32, 3, strides=1, padding='same')(x)
x = BatchNormalization()(x)
x = LeakyReLU()(x)

# 최종 output
decoder_output = Conv2DTranspose(
    1, 3, strides=1, padding='same', activation='sigmoid')(x)


decoder = Model(decoder_input, decoder_output)
decoder.summary()
###############################################################

LEARNING_RATE = 0.0005
BATCH_SIZE = 32
encoder_in = Input(shape=(100, 100, 1))
x = encoder(encoder_in)
decoder_out = decoder(x)
auto_encoder = Model(encoder_in, decoder_out)

tf.keras.utils.plot_model(encoder,to_file='encoder.png')
tf.keras.utils.plot_model(decoder,to_file='decoder.png')
tf.keras.utils.plot_model(auto_encoder,to_file='autoencoder.png')

auto_encoder.compile(optimizer=tf.keras.optimizers.Adam(
    LEARNING_RATE), loss=tf.keras.losses.MeanSquaredError())
auto_encoder.fit(x_train, y_train,
                 batch_size=BATCH_SIZE,
                 epochs=epoch
                 )
auto_encoder.save( str(epoch)+'_epochs.h5')
########################### plot ###########################
xy = encoder.predict(x_train)
validation = encoder.predict(validation_set)
plt.figure(figsize=(15, 12))
plt.scatter(x=xy[:, 0], y=xy[:, 1], c=y_train,
            cmap=plt.get_cmap('Paired'), s=3)
a=np.mean(xy[:,0])
b=np.mean(xy[:,1])
# plt.scatter(x=validation[:,0],y=validation[:,1],c=['r','b','g','m'])
plt.colorbar()
plt.show()
###############################################################