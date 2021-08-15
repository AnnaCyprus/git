import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

with open('many_genes_correlation.txt', 'r') as f:
    a = f.read().splitlines()
s = a[0]
x = s.split()
x = [float(s) for s in x]
s = a[1]
y = s.split()
y = [float(s) for s in y]

a = []
b = []
for i in range(0, 106):
    if x[i] != 0 and y[i] != 0:
        a.append(x[i])
        b.append(y[i])

a = pd.Series(a)
b = pd.Series(b)
print(a.corr(b))
