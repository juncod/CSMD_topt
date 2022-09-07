import os
from this import d
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
    X_compl.append(sample)
    X_compl.append(fliped_sample)
    for j in range(3):
        sample = np.flip(sample.transpose(), axis=1)
        fliped_sample = np.flip(fliped_sample.transpose(), axis=1)
        X_compl.append(sample)
        X_compl.append(fliped_sample)

arr = np.stack(X_compl, axis=0)

print(arr.shape)
