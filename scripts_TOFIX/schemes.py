import numpy as np
import pandas as pd

a = pd.read_csv("schemes", header=None)
a = a.drop([8], axis=1)
print(a)

b = pd.DataFrame(columns=list(range(8)))

for index, row in a.iterrows():
    mean = np.mean(row)
    b.loc[index] = row / np.mean(row)

print(b)
