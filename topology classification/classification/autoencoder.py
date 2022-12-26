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
os.environ["CUDA_VISIBLE_DEVICES"] = "2"








epoch = 1000
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

history = conv_ae.fit(x_train, x_train,batch_size=BATCH_SIZE,epochs=epoch,validation_data=(x_valid, x_valid))








X_train_compressed = conv_encoder.predict(x_train)
nsamples, nx, ny, nkernel = X_train_compressed.shape
tsne = TSNE()
X_train_2D = tsne.fit_transform(X_train_compressed.reshape((nsamples,nx*ny*nkernel)))
X_train_2D = (X_train_2D - X_train_2D.min()) / (X_train_2D.max() - X_train_2D.min())

def plot_feature(X_2D,Y):
    plt.scatter(X_2D[:, 0], X_2D[:, 1], c=Y, cmap=plt.get_cmap('Paired'), marker=",", s=0.1)
    plt.savefig('plot_feature.png',dpi=500)
    plt.show()
plot_feature(X_train_2D,y_train)
# plt.scatter(X_train_2D[:, 0], X_train_2D[:, 1], c=y_train, cmap=plt.get_cmap('Paired'), s=3)
# plt.axis("off")
# plt.show()








def plot_loss(history):
    plt.plot(history.history['loss'], label='Train loss')
    plt.plot(history.history['val_loss'], label='Validation loss')
    plt.title('Loss')
    plt.xlabel('epoch')
    plt.ylabel('loss')
    plt.legend()
    plt.ylim(0,1)
    plt.savefig('plot_loss.png',dpi=500)
    plt.show()
def plot_accuracy(history):
    plt.plot(history.history['rounded_accuracy'], label='Train accuracy')
    plt.plot(history.history['val_rounded_accuracy'], label='Validation accuracy')
    plt.title('Accuracy')
    plt.xlabel('epoch')
    plt.ylabel('accuracy')
    plt.legend()
    plt.ylim(0,1)
    plt.savefig('plot_accuracy.png',dpi=500)
    plt.show()
plot_loss(history)
plot_accuracy(history)








grid_list = []
for gl in range(0,102):
    grid_list.append(gl/100)
def plot_feature(X_2D,Y):
    plt.scatter(X_2D[:, 0], X_2D[:, 1], c=Y, cmap=plt.get_cmap('Paired'), marker=",", s=0.05)
    plt.xticks(grid_list)
    plt.yticks(grid_list)
    plt.grid(linewidth=0.5)
    plt.savefig('plot_feature.png',dpi=500)
    plt.show()
    
plot_feature(X_train_2D,y_train)







picked = []
for i in range(0,len(X_train_2D)):
    x_start = 0.22
    x_end = 0.29
    y_start = 0.5
    y_end = 0.6
    if X_train_2D[i,0] > x_start and X_train_2D[i,0] < x_end and X_train_2D[i,1] > y_start and X_train_2D[i,1] < y_end:
        picked.append(i)
print(picked)







number = 1350

pixel_plot = plt.figure()
  
# plotting a plot
pixel_plot.add_axes()
  
# customizing plot
plt.title("pixel_plot")
pixel_plot = plt.imshow(
  x_train[number-1,:,:,0], cmap='Greys', interpolation='nearest')
plt.savefig(str(number)+'.png',dpi=500)
plt.colorbar(pixel_plot)