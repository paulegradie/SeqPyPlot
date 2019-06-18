
# coding: utf-8

# # Experiments to detect heteroskedasticity



import pandas as pd
import numpy as np
import seaborn as sns
from random import randint as rand
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

from sklearn.linear_model import LinearRegression
from sklearn.metrics.pairwise import euclidean_distances

from scipy.linalg import svd

from SeqPy.container.data_container import DataContainer
from seqpyplot.parsers.config_parser import config_parser
from pathlib import Path

from matplotlib import rcParams

rcParams['figure.figsize'] = (10, 10)

pd.options.mode.chained_assignment = None 


# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Helper functions

# # Detect Heteroskedasticity



config = '../examples/example_config.ini'
config_obj = config_parser(config)




config_obj




container_obj.file_pairs




# load the data container_obj
container_obj = DataContainer(config_obj)
data, ercc_data = container_obj.parse_input()
data = container_obj.normalize_file_pairs(data) # Single df of normalized data




cols = data.columns
cols




df = data[['D3_Cont', 'D3_Treat']]
df.loc[:, 'mean'] = df.mean(axis=1)




d1 = df[['D3_Cont', 'mean']]
d2 = df[['D3_Treat', 'mean']]




fig, ax = plt.subplots()
d1.plot('mean', 'D3_Cont', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='blue', alpha=0.4)
d2.plot('mean', 'D3_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red', alpha=0.4);
ax.set_title("Raw unnormalized data");




# Quick reset for this cell
d1 = df[['D1_Cont', 'mean']]
d2 = df[['D1_Treat', 'mean']]

# define regression objects
regCont = LinearRegression(fit_intercept=True)
regTreat = LinearRegression(fit_intercept=True)

# fit regression
regCont.fit(d1['D1_Cont'].values.reshape(-1, 1), d1['mean'].values.reshape(-1, 1))
regTreat.fit(d2['D1_Treat'].values.reshape(-1, 1), d1['mean'].values.reshape(-1, 1))

print(regCont.coef_, regCont.intercept_)
print(regTreat.coef_, regTreat.intercept_)

# Correct bias
d1['D1_Cont'] = d1['D1_Cont'] - regCont.intercept_
d2['D1_Treat'] = d2['D1_Treat'] - regTreat.intercept_

# Plot regression lines
fig, ax = plt.subplots()
d1.plot('mean', 'D1_Cont', kind='scatter', ax=ax, color='blue', alpha=0.4)
d2.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 8000), ylim=(0, 8000), ax=ax, color='red',  alpha=0.4)
plt.plot([0, 8000], [0.0, regCont.coef_ * 8000], linestyle='--', color='black')
plt.plot([0, 8000], [0.0, regTreat.coef_ * 8000], linestyle='--', color='black');
ax.set_title("bias corrected, with best fit lines");

plt.ylim((-500, 8000))
plt.xlim((-500, 8000));
plt.axvline(0, linestyle='--')
plt.axhline(0, linestyle='--')
plt.legend();




plt.scatter(d1['D1_Cont'], d2['D1_Treat'])




fig, ax = plt.subplots()
x1 = d1['mean']
y1 = d1['D1_Cont'].apply(lambda x: np.log2(x))
x2 = d2['mean']
y2 = d2['D1_Treat'].apply(lambda x: np.log2(x))

ax.scatter(x1, y1, c='red', alpha=0.3)
ax.scatter(x2, y2, c='blue', alpha=0.3)
ax.set_title("Raw unnormalized data");
ax.set_xlim((0, 2500))




container_obj.file_pairs




# results = list()
# for pair in container_obj.file_pairs:
#     df = data[list(pair)]
#     mean = df.mean(axis=1)
#     df = df.sub(mean)
#     results.append(df)


# centeral all control observations on zero



results = list()
for pair in container_obj.file_pairs:
    df = data[list(pair)]
    df = df.apply(lambda x: x - df[df.columns[0]])
    results.append(df)
dat = pd.concat(results, axis=1)
dat = dat[[x[1] for x in container_obj.file_pairs]]




dat.loc[:, 'std'] = dat.std(axis=1)




dat.head()




dat.loc['Mir6957']




with open('dataout.txt', 'w+') as outfile:
    for i in top_std_genes:
        outfile.write(str(i) + '\n')




dat.sort_values(by=['std'], ascending=False)['std'].std()




top_std_genes = dat[dat.sort_values(by=['std'], ascending=False)['std'] >  dat.sort_values(by=['std'], ascending=False)['std'].std()].sort_values(by=['std'], ascending=False).index




len(top_std_genes)




dat['std'].plot(kind='hist')




np.random.randint(0, 100, 20).std()




dat




dat.head(20).std(axis=1
                 
                )


