import tensorflow as tf
from tensorflow import keras 
from tensorflow.keras.models import Model, Sequential, load_model
from tensorflow.keras.layers import Conv2D, Conv2DTranspose, SeparableConv2D, BatchNormalization, Reshape, UpSampling2D, MaxPool2D,SeparableConv2D
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE 
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
import umap
import os

epoch = 1000
LEARNING_RATE = 0.001
BATCH_SIZE = 1024

compl_path = './compl'
compl_file_list = os.listdir(compl_path)
compl_file_list_py = [
    file for file in compl_file_list if file.endswith('.xlsx')]

heat_path = './heat'
heat_file_list = os.listdir(heat_path)
heat_file_list_py = [file for file in heat_file_list if file.endswith('.xlsx')]


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
    


X_list = X_compl + X_heat
x_train = np.stack(X_list, axis=0)
x_train = x_train.reshape(-1, 128, 128, 1)

y_train = []
Y_compl = ['b' for i in range(len(X_compl))]
Y_heat = ['r' for i in range(len(X_heat))]
y_train = Y_compl + Y_heat
y_train = np.array(y_train)

train_noise = np.random.normal(loc=0.0, scale=0.5, size=x_train.shape)
x_train_noisy = np.clip(x_train + train_noise, 0, 1)


def MSE(y_true, y_pred):
    return keras.metrics.MeanSquaredError(y_true, y_pred)
def PSNR(y_true, y_pred):
    return tf.image.psnr(y_true, y_pred, max_val=1.0)
def SSIM(y_true, y_pred):
    return tf.image.ssim(y_true, y_pred, max_val=1.0)
conv_encoder = load_model('encoder.h5',custom_objects={"PSNR": PSNR ,"SSIM":SSIM, "MSE":MSE})

conv_encoder.compile(loss=tf.keras.losses.BinaryCrossentropy(), optimizer=keras.optimizers.Adam(LEARNING_RATE),
                metrics=[tf.keras.metrics.MeanSquaredError(),PSNR,SSIM])


X_train_compressed = conv_encoder.predict(x_train)
nsamples, nx, ny, nkernel = X_train_compressed.shape
reducer = umap.UMAP(n_components=2)
X_train_2D = reducer.fit_transform(X_train_compressed.reshape((nsamples,nx*ny*nkernel)))
X_train_2D = (X_train_2D - X_train_2D.min()) / (X_train_2D.max() - X_train_2D.min())
y_train = []
Y_compl = ['b' for i in range(len(X_compl))]
Y_heat = ['r' for i in range(len(X_heat))]
y_train = Y_compl + Y_heat
y_train = np.array(y_train)
def plot_feature(X_2D,Y):
    plt.title('TSNE_noise')
    plt.scatter(X_2D[:, 0], X_2D[:, 1], c=Y, marker=",", s=1)
    plt.savefig('TSNE.png',dpi=500)
    plt.show()
plot_feature(X_train_2D,y_train)