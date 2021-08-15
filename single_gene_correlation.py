import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

with open('single_gene_correlation.txt', 'r') as f:
    a = f.read().splitlines()
s = a[len(a) - 1]
y_lab = s.split()
a.pop()
z_values = []
for s in a:
    b = s.split()
    b = [float(s1) for s1 in b]
    z_values.append(b)

a = []
for i in range(0, 106):
    n = 0
    for j in range(0, len(z_values)):
        if z_values[j][i] != 0:
            n = 1
    if n != 0:
        b = []
        for j in range(0, len(z_values)):
            b.append(z_values[j][i])
        a.append(b)
    
df = pd.DataFrame(a, columns = y_lab)
corr_mat = df.corr()
plt.figure(figsize=(10, 7))
sns.set(font_scale = 1)
plt.title('Correlation')
sns.heatmap(corr_mat, annot=True, cmap='vlag', vmin=-1, vmax=1)
plt.xticks(rotation = 0)
plt.yticks(rotation = 0)
plt.show()
