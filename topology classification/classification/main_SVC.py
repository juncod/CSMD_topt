import os
import pandas as pd
from sklearn.svm import SVC
import numpy as np

compl_path = './compl'
compl_file_list = os.listdir(compl_path)
compl_file_list_py = [
    file for file in compl_file_list if file.endswith('.xlsx')]

heat_path = './heat'
heat_file_list = os.listdir(heat_path)
heat_file_list_py = [file for file in heat_file_list if file.endswith('.xlsx')]

test_path = './test'
test_file_list = os.listdir(test_path)
test_file_list_py = [file for file in test_file_list if file.endswith('.xlsx')]
print(test_file_list_py)

X_compl = []
for i in compl_file_list_py:
    filename = compl_path+'/'+i
    df = pd.read_excel(filename, engine='openpyxl',
                       header=None, index_col=None)
    sample = df.to_numpy()
    fliped_sample = np.flip(df.to_numpy(), axis=1)
    X_compl.append(sample.ravel())
    X_compl.append(fliped_sample.ravel())
    for j in range(3):
        sample = np.flip(sample.transpose(), axis=1)
        fliped_sample = np.flip(fliped_sample.transpose(), axis=1)
        X_compl.append(sample.ravel())
        X_compl.append(fliped_sample.ravel())

X_heat = []
for i in heat_file_list_py:
    filename = heat_path+'/'+i
    df = pd.read_excel(filename, engine='openpyxl',
                       header=None, index_col=None)
    sample = df.to_numpy()
    fliped_sample = np.flip(df.to_numpy(), axis=1)
    X_heat.append(sample.ravel())
    X_heat.append(fliped_sample.ravel())
    for j in range(3):
        sample = np.flip(sample.transpose(), axis=1)
        fliped_sample = np.flip(fliped_sample.transpose(), axis=1)
        X_heat.append(sample.ravel())
        X_heat.append(fliped_sample.ravel())

X_test = []
for i in test_file_list_py:
    filename = test_path+'/'+i
    df = pd.read_excel(filename, engine='openpyxl',
                       header=None, index_col=None)
    sample = df.to_numpy()
    X_test.append(sample.ravel())
print(len(X_compl))

Y_compl = [0 for i in range(len(X_compl))]
Y_heat = [1 for i in range(len(X_heat))]

X_sample = X_compl + X_heat
Y_sample = Y_compl + Y_heat
classifier = SVC(kernel='sigmoid')  # linear poly rbf sigmoid
classifier.fit(X_sample, Y_sample)


for i in range(len(X_test)):
    print(classifier.predict(X_test[i-1].reshape(1, -1)))
