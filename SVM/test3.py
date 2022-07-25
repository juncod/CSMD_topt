import numpy as np

arr_2d = np.arange(25).reshape(5, 5)
b = np.flip(arr_2d.transpose(), axis=1)
c = np.flip(b.transpose(), axis=1)
print(arr_2d, '\n', b, '\n', c)
