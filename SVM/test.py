import os
import pandas as pd
from sklearn.svm import SVC

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

X_sample = []
for i in compl_file_list_py:
    filename = compl_path+'/'+i
    df = pd.read_excel(filename, engine='openpyxl',
                       header=None, index_col=None)
    sample = df.to_numpy()
    X_sample.append(sample.ravel())
    X_sample.append(sample[::-1].ravel())

    break
