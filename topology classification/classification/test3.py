import tensorflow as tf
from tensorflow import keras 
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Conv2D, Conv2DTranspose, Dense, Flatten, Dropout, BatchNormalization, Reshape, LeakyReLU, MaxPool2D
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
from tensorflow.keras.models import load_model
from sklearn.manifold import TSNE
import os
# os.environ["CUDA_VISIBLE_DEVICES"] = "2"

epoch = 1
LEARNING_RATE = 0.0002
BATCH_SIZE = 1024

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

valid_list = []
for i in valid_file_list_py:
    filename = valid_path+'/'+i
    df = pd.read_excel(filename, engine='openpyxl',
                       header=None, index_col=None)
    sample = df.to_numpy()
    valid_list.append(sample)

x_valid = np.stack(valid_list, axis=0)
x_valid = x_valid.reshape(-1, 100, 100, 1)
    


X_list = X_compl + X_heat
x_train = np.stack(X_list, axis=0)
x_train = x_train.reshape(-1, 100, 100, 1)

y_train = []
Y_compl = [0 for i in range(len(X_compl))]
Y_heat = [1 for i in range(len(X_heat))]
y_train = Y_compl + Y_heat
y_train = np.array(y_train)





def rounded_accuracy(y_true, y_pred):
    return keras.metrics.binary_accuracy(tf.round(y_true), tf.round(y_pred))
####################### Model ###################################
conv_encoder = Sequential([
    Reshape([100,100,1],input_shape=[100,100]),
    Conv2D(16,kernel_size=5,padding="same",activation="selu"),
    Conv2D(16,kernel_size=5,padding="same",activation="selu"),
    MaxPool2D(pool_size=5),
    Conv2D(32,kernel_size=3,padding="same",activation="selu"),
    Conv2D(32,kernel_size=3,padding="same",activation="selu"),
    MaxPool2D(pool_size=2)
])
conv_decoder = Sequential([
    Conv2DTranspose(16,kernel_size=3,strides=1,padding="same",activation="selu"),
    Conv2DTranspose(16,kernel_size=3,strides=2,padding="same",activation="selu"),
    Conv2DTranspose(1,kernel_size=5,strides=5,padding="same",activation="selu"),
    Conv2DTranspose(1,kernel_size=5,strides=1,padding="same",activation="sigmoid"),
    Reshape([100,100])
])
conv_ae = Sequential([conv_encoder,conv_decoder])

conv_ae.compile(loss="binary_crossentropy", optimizer=keras.optimizers.SGD(LEARNING_RATE),
                metrics=[rounded_accuracy])
# auto_encoder.compile(optimizer=tf.keras.optimizers.Adam(
#     LEARNING_RATE), loss=tf.keras.losses.MeanSquaredError())

###############################################################
conv_encoder.summary()
conv_decoder.summary()
conv_ae.summary()

conv_ae.fit(x_train, x_train,batch_size=BATCH_SIZE,epochs=epoch,validation_data=(x_valid, x_valid)) # batch_size=BATCH_SIZE,
# conv_ae.save( str(epoch)+'_epochs.h5')


X_train_compressed = conv_encoder.predict(x_train)
nsamples, nx, ny, nkernel = X_train_compressed.shape
tsne = TSNE()
X_train_2D = tsne.fit_transform(X_train_compressed.reshape((nsamples,nx*ny*nkernel)))
X_train_2D = (X_train_2D - X_train_2D.min()) / (X_train_2D.max() - X_train_2D.min())

plt.scatter(X_train_2D[:, 0], X_train_2D[:, 1], c=y_train, cmap=plt.get_cmap('Paired'), s=3)
plt.axis("off")
plt.show()





########################### plot ###########################
# xy = conv_encoder.predict(x_train)
# validation = conv_encoder.predict(validation_set)
# plt.figure(figsize=(15, 12))
# plt.scatter(x=xy[:, 0], y=xy[:, 1], c=y_train,
#             cmap=plt.get_cmap('Paired'), s=3)

# plt.scatter(x=validation[:,0],y=validation[:,1],c=['r','b','g','m'])
# plt.colorbar()
# plt.show()
###############################################################

# plt.figure(figsize=(10,8))
# cmap = plt.cm.tab10
# plt.scatter(X_valid_2D[:, 0], X_valid_2D[:, 1], c=y_train, s=10, cmap=cmap)
# image_positions = np.array([[1., 1.]])
# for index, position in enumerate(X_valid_2D):
#     dist = np.sum((position - image_positions) ** 2, axis=1)
#     if np.min(dist) > 0.02: # if far enough from other images
#         image_positions = np.r_[image_positions, [position]]
#         imagebox = mpl.offsetbox.AnnotationBbox(
#             mpl.offsetbox.OffsetImage(x_train[index], cmap="binary"),
#             position, bboxprops={"edgecolor": cmap(y_train[index]), "lw": 1})
#         plt.gca().add_artist(imagebox)
# plt.axis("off")
# plt.show()