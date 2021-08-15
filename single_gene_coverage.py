import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

x_lab = []
for i in range(1, 107):
    x_lab.append(i)
with open('single_gene_coverage.txt', 'r') as f:
    a = f.read().splitlines()
s = a[len(a) - 1]
y_lab = s.split()
a.pop()
z_values = []
for s in a:
    b = s.split()
    b = [float(s1) for s1 in b]
    z_values.append(b)

df = pd.DataFrame(z_values, columns = x_lab, index = y_lab)
plt.figure(figsize=(10, 7))
sns.set(font_scale = 0.5)
CDRs = [27, 38, 56, 65, 105]
for i in range(0, len(CDRs)):
    if i % 2 == 0:
        xc = CDRs[i] - 1
    else:
        xc = CDRs[i]
    plt.axvline(x = xc, color = 'r')
fam = [4, 10, 11, 9, 2]
yc = 0
for i in range(0, len(fam) - 1):
    yc += fam[i]
    plt.axhline(y = yc, color = 'r')
plt.title('Coverage')
sns.heatmap(df, annot = False, cmap="Blues", vmin=0, vmax=1)
plt.show()
