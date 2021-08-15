import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

x = []
for i in range(1, 107):
    x.append(i)
x = np.array(x)
with open('many_genes_plot.txt', 'r') as f:
    a = f.read().splitlines()
s = a[0]
y = s.split()
y = [float(s) for s in y]
y1 = np.array(y)
s = a[1]
y = s.split()
y = [float(s) for s in y]
y2 = np.array(y)

plt.figure(figsize=(10, 7))
plt.xlim([1, 106])
plt.ylim([0, 1])
plt.plot(x, y1, 'b', label = 'germline')
plt.plot(x, y2, 'r', label = 'SHM')
CDRs = [27, 38, 56, 65, 105]
for xc in CDRs:
    plt.axvline(x = xc, color = 'g')
plt.xlabel('IMGT aminoacid position')
plt.title('Allelic Variation and Somatic Hypermutation')
legend = plt.legend(loc = 'upper right')
plt.show()
