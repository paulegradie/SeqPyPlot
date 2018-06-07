
# coding: utf-8

# # Heteroskedasticity and Bias correction using linear transformations
# 
# This notebook works through an idea for correcting heteroskedasticity which uses rotation matrices to align samples. Previous iterations of seqpyplot (SPPLOT v0.3 and below) output a scatter plot that plots the expression values of control and treated samples and against their means. These plots revealed what appeared to be a bias where one sample tended towards overall higher expression than the other. Given the type of data being analyzed (Time series RNA seq data), it seemedunlikely that there would be such a bias inherent to the experiment.
# 
# 
# The steps to performing this transformation are as follows
# 1. For a given sample pair, compute the mean. 
# 2. For each sample and its repestive mean computed in step 1, perform linear regression to compute a regression coefficient (a line of best fit) through the sample/mean pair as well as a bias.
# 3. Zero out the sample bias by addition or subtraction across all points.
# 4. If the regression coeffieicients do not match, use the coefficients to compute the angle difference.
# 5. Compute a rotation matrix to rotate the data about the origin until they match, resulting in proper alignment between the control and treated samples.
# 

# References:
#     Rotation Matrices: https://en.wikipedia.org/wiki/Rotation_matrix



import pandas as pd
import numpy as np
from random import randint as rand
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

from sklearn.linear_model import LinearRegression
from sklearn.metrics.pairwise import euclidean_distances

from scipy.linalg import svd

from seqpyplot.container.data_container import DataContainer
from seqpyplot.parsers.config_parser import config_parser
from pathlib import Path

from matplotlib import rcParams

rcParams['figure.figsize'] = (15, 15)

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



def calc_theta(coef1, coef2):
    "Returns an angle in radians"
    return np.abs(
       np.arctan(np.abs(coef1 - coef2) / (1. + (coef1 * coef2)))
   )




def compute_rot_mat(rad, coef=.5):
    " Compute a rotation matrix using rad for a given regression coefficient "
    if coef < 1.0:
        rotation_matrix = np.array([[np.cos(rad), -np.sin(rad)],
                                    [np.sin(rad), np.cos(rad)]])
    else:
        rotation_matrix = np.array([[np.cos(rad), np.sin(rad)],
                                    [-np.sin(rad), np.cos(rad)]])  
    return rotation_matrix


# # Simple example of line rotation
# Example of rotating a line around the origin to match the slope of another line



slope1 = 1.1
slope2 = 2.0

line1 = np.array([slope1 * x for x in range(10)])
line2 = np.array([slope2 * x for x in range(10)])

xs = list(range(10))


# #### Original Lines



# Plot lines
plt.plot(xs, line1, color='black', label='line to rotate');
plt.plot(xs, line2, color='red', linewidth=5, label='stationary');
plt.axvline(0, linestyle='--')
plt.axhline(0, linestyle='--');

plt.annotate(s=('<-- we\'ll rotate this line'), xy=(xs[3]+0.08, line1[3]), rotation=-36)
plt.annotate(s=('<-- Direction of rotation'), xy=(xs[3] + 0.75, line1[3]+3.5), rotation=-15)

plt.ylim((-1, 10))
plt.xlim((-1, 10))
plt.legend(loc='upper right');


# #### Rotated lines



# Compute angle
angle_diff = calc_theta(slope1, slope2)
angle_diff




# Compute rotation matrix
rot_matrix = compute_rot_mat(angle_diff)
rot_matrix




# rotate line 1 (black line)
new_line1 = list()
for x, y in zip(xs, line1):
    # need shape [[#], [#]]
    old_point = np.array([[x], [y]])
    new_point = np.dot(rot_matrix, old_point)
    new_line1.append(new_point)
new_line1 = np.squeeze(np.asarray(new_line1))




xs[6], line1[6]




plt.plot(xs, line2, color='red', linewidth=5, alpha=0.7, label='Stationary');

plt.plot(xs, line1, color='black', label='line to rotate')
plt.scatter(new_line1[:, 0], new_line1[:, 1], color='black', s=95)
plt.plot(new_line1[:, 0], new_line1[:, 1], color='black', linestyle='--', label='rotated line')

plt.annotate(s=('Original Line'), xy=(xs[6] + 0.7, line1[6]+0.3))
plt.annotate(s=('<-- Direction of rotation'), xy=(xs[6] - 1.7, line1[6]+0.8), rotation=-36)

plt.axvline(0, linestyle='--')
plt.axhline(0, linestyle='--')

for (x1, y1), (x2, y2) in zip(zip(xs, line1), new_line1):
    plt.plot([x1, x2], [y1, y2], linestyle='--', color='black', alpha=0.4)

plt.ylim((-1, 10))
plt.xlim((-1, 10));
plt.legend(loc='upper right');


# # Example using raw expression Data
# (No TMM normalization)



config = '../examples/example_config.ini'
config_obj = config_parser(config)




# load the data container_obj
container_obj = DataContainer(config_obj)
data, ercc_data = container_obj.parse_input()




cols = data.columns
cols




df = data[['D1_Cont', 'D1_Treat']]
df.loc[:, 'mean'] = df.mean(axis=1)




df.head()




d1 = df[['D1_Cont', 'mean']]
d2 = df[['D1_Treat', 'mean']]




fig, ax = plt.subplots()
d1.plot('mean', 'D1_Cont', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='blue', alpha=0.4)
d2.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red', alpha=0.4);
plt.annotate('The bias between samples is clearly seen in this plot!', (500, 3900), fontsize=14)
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


# #### Compute rotation linear transformation



correction_theta = calc_theta(np.squeeze(regCont.coef_), np.squeeze(regTreat.coef_))
correction_theta  # in radians

rotation_matrix = compute_rot_mat(correction_theta)
rotation_matrix

new_treat = np.array([np.dot(rotation_matrix, d2.values[i, :]) for i in range(len(d2.values))])
new_treat

d2_cor = d2.copy()
d2_cor.loc[:, 'D1_Treat'] = new_treat[:, 0]
d2_cor.loc[:, 'mean'] = new_treat[:, 1]

fig, ax = plt.subplots()
d1.plot('mean', 'D1_Cont', kind='scatter', xlim=(0, 20000), ylim=(0, 20000), ax=ax, color='blue', alpha=0.4)
d2_cor.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red', alpha=0.4);
ax.set_title("No TMM, linearly transformed")

plt.ylim((-1000, 8000))
plt.xlim((-1000, 8000));
plt.axvline(0, linestyle='--')
plt.axhline(0, linestyle='--')
plt.legend();


# Successs! We have eliminated the sample divergence!

# 
# # Use TMM before linear transformation
# 
# There could still be advantage in normalizing the samples prior to performing the linear transformation.



# load the data container_obj
config = '../examples/example_config.ini'
config_obj = config_parser(config)

container_obj = DataContainer(config_obj)
data, ercc_data = container_obj.parse_input()

data = container_obj.normalize_file_pairs(data) # Single df of normalized data

normed_df = data[['D1_Cont', 'D1_Treat']].copy()
normed_df.loc[:, 'mean'] = normed_df.mean(axis=1)
 
regCont = LinearRegression(fit_intercept=True)
regTreat = LinearRegression(fit_intercept=True)

regCont.fit(normed_df['D1_Cont'].values.reshape(-1, 1), normed_df['mean'].values.reshape(-1, 1))
regTreat.fit(normed_df['D1_Treat'].values.reshape(-1, 1), normed_df['mean'].values.reshape(-1, 1))

normed_df['D1_Cont'] = normed_df['D1_Cont'] - regCont.intercept_
normed_df['D1_Treat'] = normed_df['D1_Treat'] - regTreat.intercept_

fig, ax = plt.subplots()
normed_df.plot('mean', 'D1_Cont', kind='scatter',
               xlim=(0, 5000), ylim=(0, 8000), ax=ax, color='blue', alpha=0.4, label='Control')
normed_df.plot('mean', 'D1_Treat', kind='scatter',
               xlim=(0, 5000), ylim=(0, 8000), ax=ax, color='red', alpha=0.4, label='Treated')

# plot regression lines, with color switch!
plt.plot([0, 8000], [0.0, regCont.coef_ * 8000], linestyle='--', color='black')
plt.plot([0, 4000], [0.0, regTreat.coef_ * 4000], linestyle='--', color='white');
plt.plot([4000, 8000], [regTreat.coef_[0] * 4000.0, regTreat.coef_[0] * 8000], linestyle='--', color='black');

ax.set_title("TMM normalized expression data");

plt.ylim((-500, 8000))
plt.xlim((-500, 8000));
plt.axvline(0, linestyle='--')
plt.axhline(0, linestyle='--')
plt.legend();


# #### Compute linear transformation



correction_theta = calc_theta(np.squeeze(regCont.coef_), np.squeeze(regTreat.coef_))
correction_theta  # in radians

rotation_matrix = compute_rot_mat(correction_theta, regTreat.coef_)
rotation_matrix

new_treat = np.array([np.dot(rotation_matrix, normed_df[['D1_Treat', 'mean']].values[i, :]) for i in range(len(normed_df))])

corr_df = normed_df.copy()
corr_df.loc[:, 'D1_Treat'] = new_treat[:, 0]
# corr_df.loc[:, 'mean'] = normed_df['mean'].values
corr_df.loc[:, 'mean'] = new_treat[:, 1]

fig, ax = plt.subplots()
normed_df.plot('mean', 'D1_Cont', kind='scatter', xlim=(0, 20000), ylim=(0, 20000), ax=ax, color='blue', alpha=0.4)
corr_df.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red', alpha=0.4, s=10);
ax.set_title("With TMM, linearly transformed");

plt.ylim((-500, 15000))
plt.xlim((-500, 15000));
plt.axvline(0, linestyle='--')
plt.axhline(0, linestyle='--')
plt.legend();




correction_theta = calc_theta(np.squeeze(regCont.coef_), np.squeeze(regTreat.coef_))
correction_theta  # in radians

rotation_matrix = compute_rot_mat(correction_theta, regTreat.coef_)
rotation_matrix

new_treat = np.array([np.dot(rotation_matrix, normed_df[['D1_Treat', 'mean']].values[i, :]) for i in range(len(normed_df))])

corr_df = normed_df.copy()
corr_df.loc[:, 'D1_Treat'] = new_treat[:, 0]
corr_df.loc[:, 'mean'] = new_treat[:, 1]




fig, ax = plt.subplots()
normed_df.plot('mean', 'D1_Cont', kind='scatter', xlim=(0, 20000), ylim=(0, 20000), ax=ax, color='blue', alpha=0.4)
corr_df.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red', alpha=0.4, s=10);
ax.set_title("With TMM, linearly transformed");

plt.ylim((-50, 500))
plt.xlim((-50, 500));
plt.axvline(0, linestyle='--')
plt.axhline(0, linestyle='--')
plt.legend();




corr_df.loc[:, 'mean'] = normed_df['mean'].values


# # Outlier detection
# 
# With the transformed data, we need a new metric for detecting outliers. While previously we described points with a single unique value and a common shared value, this transformed data is now represented with four unique values (two points for each sample).
# 
# Fortunately, we can still measure the distance between these points using euclidean distance. We can also create a new regression line using all of this new data from which we can establish a logfold threshold - thus providing a baseline for euclidean distance measurements. We can compute this threshold distance using the following algorithm:
# 
# 1. Compute the regression coefficient for the new data (concat the data frames for each sample).
# 2. Compute the upper and lower thersholds given the user selected logfold threshold.
# 3. Compute the euclidean distance between all sample point pairs.



normed_df.drop('D1_Treat', axis=1, inplace=True)
normed_df.rename(columns={'mean': 'cont_mean'}, inplace=True)
normed_df.head()




corr_df.drop('D1_Cont', axis=1, inplace=True)
corr_df.rename(columns={'mean': 'treat_mean'}, inplace=True)
corr_df.head()




def euclidean_distance(p, q):
    return np.sqrt(
        np.sum(
    ((p[0] - q[0]) ** 2) + ((p[1] - q[1]) ** 2)
    ))




distances = list()
for i in normed_df.index:
    p = normed_df.loc[i, :].values
    q = corr_df.loc[i, :].values
    distances.append(euclidean_distance(p, q))




filt_points = sorted(list(filter(lambda x: (x < 1000) & (x > 50), distances)))
plt.hist(filt_points[50:], bins=50);


# #### Plot a log2Fold range around the regression line



def calc_bounds(series_mean, var=0.5):
    """
    Use some linear algebra to calculate log2fold range around a given value based on user logfold parameter
    """
    return np.array([i - ((2.0 * i) / ((2.0 ** var) + 1)) for i in series_mean])




main = np.array([0, 5, 10])




spread = calc_bounds(main)




spread




plt.plot(main)
plt.plot(main + spread)
plt.plot(main - spread);


# Plot Values with min threshhold



low = 1984.0
hi = 1984.2

tot_data = pd.concat([normed_df.rename(columns={'D1_Cont': 'expr', 'cont_mean': 'mean'}, inplace=False),
                      corr_df.rename(columns={'D1_Treat': 'expr', 'treat_mean': 'mean'}, inplace=False)], axis=0)

tot_reg = LinearRegression(fit_intercept=True)
tot_reg.fit(tot_data['mean'].values.reshape(-1, 1), tot_data['expr'].values.reshape(-1, 1))

series_mean = np.array([x * tot_reg.coef_[0] for x in range(int(hi))])

spread = calc_bounds(series_mean, var=0.4)

plt.scatter(
    np.abs(normed_df[(normed_df.values[:, 1] > low) & (normed_df.values[:, 1] < hi)].values[:, 1]),
    np.abs(normed_df[(normed_df.values[:, 1] > low) & (normed_df.values[:, 1] < hi)].values[:, 0]),
    color='blue',
    s=10,
    alpha=0.3)

plt.scatter(
    np.abs(corr_df[(corr_df.values[:, 1] > low) & (corr_df.values[:, 1] < hi)].values[:, 1]),
    np.abs(corr_df[(corr_df.values[:, 1] > low) & (corr_df.values[:, 1] < hi)].values[:, 0]),
    color='red',
    s=10,
    alpha=0.3)
bias = 5

plt.plot(np.abs(series_mean + spread + bias), color='black', linestyle='--')
plt.plot(np.abs(series_mean - spread - bias), color='black', linestyle='--')

# plt.ylim((low - 50, hi + 1000))
plt.xlim((low - 500, hi + 500));
plt.axvline(0, linestyle='--')
plt.axhline(0, linestyle='--')

plt.plot([0, hi], [0, tot_reg.coef_[0][0] * hi], color='black', linewidth=3, linestyle='--')
plt.legend();




normed_df[(normed_df.values[:, 1] > low) & (normed_df.values[:, 1] < hi)]




corr_df[(corr_df.values[:, 1] > low) & (corr_df.values[:, 1] < hi)]




normed_df.loc['Tead1']




corr_df.loc['Tead1']




(data.loc['Tead1'][0] + data.loc['Tead1'][3]) / 2.




low = 1983
hi = 1985

tot_data = pd.concat([normed_df.rename(columns={'D1_Cont': 'expr', 'cont_mean': 'mean'}, inplace=False),
                      corr_df.rename(columns={'D1_Treat': 'expr', 'treat_mean': 'mean'}, inplace=False)], axis=0)

tot_reg = LinearRegression(fit_intercept=True)
tot_reg.fit(tot_data['mean'].values.reshape(-1, 1), tot_data['expr'].values.reshape(-1, 1))

series_mean = np.array([x * tot_reg.coef_[0] for x in range(hi)])

spread = calc_bounds(series_mean, var=0.4)

plt.scatter(
normed_df.loc['Tead1'].values[0], normed_df.loc['Tead1'].values[1]
)

plt.scatter(
corr_df.loc['Tead1'].values[0], corr_df.loc['Tead1'].values[1]
)

x = data.loc['Tead1'].values[0]
y = data.loc['Tead1'].values[3]
meen = (x + y) / 2

plt.scatter(
meen, x,
color='red',
s=90)

plt.scatter(
meen, y,
color='red',
s=90)

bias = 5

plt.plot(np.abs(series_mean + spread + bias), color='black', linestyle='--')
plt.plot(np.abs(series_mean - spread - bias), color='black', linestyle='--')

# plt.ylim((low - 50, hi + 1000))
plt.xlim((low - 500, hi + 500));
plt.axvline(0, linestyle='--')
plt.axhline(0, linestyle='--')

plt.plot([0, hi], [0, tot_reg.coef_[0][0] * hi], color='black', linewidth=3, linestyle='--')
plt.legend();




x = data.loc['Tead1'].values[0]
y = data.loc['Tead1'].values[3]
meen = (x + y) / 2




x, y, meen




(data.loc['Tead1'][0] + data.loc['Tead1'][3]) / 2.




norm_x = normed_df.loc['Tead1'][0]
corr_y = corr_df.loc['Tead1'][0]
orig_mean = normed_df.loc['Tead1'][1]




norm_x, corr_y, orig_mean




data.loc['Tead1']




low = 1983
hi = 1985

tot_data = pd.concat([normed_df.rename(columns={'D1_Cont': 'expr', 'cont_mean': 'mean'}, inplace=False),
                      corr_df.rename(columns={'D1_Treat': 'expr', 'treat_mean': 'mean'}, inplace=False)], axis=0)

tot_reg = LinearRegression(fit_intercept=True)
tot_reg.fit(tot_data['mean'].values.reshape(-1, 1), tot_data['expr'].values.reshape(-1, 1))

series_mean = np.array([x * tot_reg.coef_[0] for x in range(hi)])

spread = calc_bounds(series_mean, var=0.4)

plt.scatter(
normed_df.loc['Tead1'].values[0], normed_df.loc['Tead1'].values[1]
)

plt.scatter(
corr_df.loc['Tead1'].values[0], corr_df.loc['Tead1'].values[1]
)

norm_x = normed_df.loc['Tead1'][0]
corr_y = corr_df.loc['Tead1'][0]
orig_mean = normed_df.loc['Tead1'][1]

plt.scatter(
orig_mean, norm_x,
color='red',
s=90)

plt.scatter(
orig_mean, corr_y,
color='red',
s=90)

bias = 5

plt.plot(np.abs(series_mean + spread + bias), color='black', linestyle='--')
plt.plot(np.abs(series_mean - spread - bias), color='black', linestyle='--')

# plt.ylim((low - 50, hi + 1000))
plt.xlim((low - 500, hi + 500));
plt.axvline(0, linestyle='--')
plt.axhline(0, linestyle='--')

plt.plot([0, hi], [0, tot_reg.coef_[0][0] * hi], color='black', linewidth=3, linestyle='--')
plt.legend();




(data.loc['Tead1'][0] + data.loc['Tead1'][3]) / 2.




joined_data.head()


# #### Compute euclidean distance betwen sample points



joined_data = pd.concat([normed_df, corr_df], axis=1)

results = list()
for row in joined_data.iterrows():
    row = row[1]
    p = [row['cont_mean'], row['D1_Cont']]
    q = [row['treat_mean'], row['D1_Treat']]
    results.append(euclidean_distance(p, q))

joined_data.loc[:, 'distance'] = results




joined_data.head()




up = 1000
low = 500

# jd = joined_data[(joined_data['distance'] > up) & (joined_data['distance'] > low)].copy()
# jd = joined_data.loc[joined_data.index[4: 5]]
jd = joined_data
plt.scatter(jd.values[:, 1], jd.values[:, 0], color='blue', alpha=0.4)
plt.scatter(jd.values[:, 3], jd.values[:, 2], color='red',  alpha=0.4)

plt.plot(series_mean + spread + bias, color='black', linestyle='--')
plt.plot(series_mean - spread - bias, color='black', linestyle='--')

plt.ylim((-50, 1000))
plt.xlim((-50, 1000));
plt.axvline(0, linestyle='--')
plt.axhline(0, linestyle='--');


# #### Computer center point for each sample pair



jd = joined_data.loc[joined_data.index[4: 5]]
jd




x1, x2 = jd.values[0][1], jd.values[0][3]
y1, y2 = jd.values[0][0], jd.values[0][2]




avex = (x1 + x2) / 2.
avey = (y1 + y2) / 2.




def center_point(p, q):
    avex = (p[0] + q[0]) / 2.
    avey = (p[1] + q[1]) / 2.
    return avex, avey




avex, avey = center_point([x1, y1], [x2, y2])




plt.scatter(jd.values[:, 1], jd.values[:, 0], color='blue', alpha=0.4)
plt.scatter(jd.values[:, 3], jd.values[:, 2], color='red',  alpha=0.4)
plt.scatter(avex, avey, s=180);




len(centers)




len(joined_data)




centers = list()
for row in joined_data.iterrows():
    row = row[1]
    p = [row['cont_mean'], row['D1_Cont']]
    q = [row['treat_mean'], row['D1_Treat']]
    centers.append(center_point(p, q))




joined_data.loc[:, 'centers'] = list(centers)




x_vals = np.array([list(x) for x in joined_data['centers'].values])[:, 0]




center_dist = calc_bounds(x_vals, var=0.4)




joined_data.loc[:, 'min_dist'] = center_dist*2




joined_data.head()




test_filtered = joined_data[joined_data['min_dist'] < joined_data['distance']]
test_filtered.head()


# # Choosing thresholds for outlier cutoffs
# 
# Before performing rotation correction, its a good idea to eliminate some of the worst outlier offenders before performing linear regression. To do this, we can calculate summary statistics on the sample pair, and 



config = '../examples/example_config.ini'
config_obj = config_parser(config)

# load the data container_obj
container_obj = DataContainer(config_obj)
data, ercc_data = container_obj.parse_input()

data = container_obj.normalize_file_pairs(data) # Single df of normalized data

df = data[['D1_Cont', 'D1_Treat']].copy()
df.loc[:, 'mean'] = df.mean(axis=1)
df.loc[:, 'var'] = (df['D1_Cont'] - df['mean']).abs()
df.loc[:, 'diff'] = (df['D1_Cont'] - df['D1_Treat']).abs()




described = df[df['diff'] > 0.0].describe(percentiles=[x /10. for x in range(10)])




deef = described.copy()
deef.drop(['count', 'mean', 'std', 'min', 'max'], axis=0, inplace=True)
columns = deef.columns
idx = deef.index.tolist()




deef


# Using the 75th percentile is typical practice. Based on thse summary stats, it doens't seem unreasonable to take up to the 90th percentile, however a conservative estimate is likely to improve results.



fig = plt.figure()
grid = plt.GridSpec(2, 6, wspace=0.4, hspace=0.3)

ax1 = plt.subplot(grid[0, 0:3])
ax1.plot(idx, deef['D1_Cont'])
plt.title('Control Counts')
plt.axvline(7.5, linestyle='--', alpha=0.4)

ax2 = plt.subplot(grid[0, 3:6], sharey=ax1)
ax2.plot(idx, deef['D1_Treat'])
plt.title('Treated Counts')
plt.axvline(7.5, linestyle='--', alpha=0.4)

ax3 = plt.subplot(grid[1, 0:2])
ax3.plot(idx, deef['mean'])
plt.title('Mean')
plt.axvline(7.5, linestyle='--', alpha=0.4)

ax4 = plt.subplot(grid[1, 2:4], sharey=ax3)
ax4.plot(idx, deef['var'])
plt.title('Var')
plt.axvline(7.5, linestyle='--', alpha=0.4)

ax5 = plt.subplot(grid[1, 4:6], sharey=ax3)
ax5.plot(idx, deef['diff'])
plt.title('Diff')
plt.axvline(7.5, linestyle='--', alpha=0.4)

plt.tight_layout();




perc_filtered = df #[df['diff'] < deef.describe().loc['75%']['diff']]




def rotate(rotation_matrix, points):
    return np.array([np.dot(rotation_matrix, points.values[i, :]) for i in range(len(points.values))])

def calc_theta(coef1, coef2):
    "returns radians"
    return np.abs(
       np.arctan(np.abs(coef1 - coef2) / (1. + (coef1 * coef2)))
   )

def compute_rot_mat(rad, coef):
    if coef < 1.0:
        rotation_matrix = np.array([[np.cos(rad), -np.sin(rad)],
                                    [np.sin(rad), np.cos(rad)]])
    else:
        rotation_matrix = np.array([[np.cos(rad), np.sin(rad)],
                                    [-np.sin(rad), np.cos(rad)]])  
    return rotation_matrix




# Compute regression
regCont = LinearRegression(fit_intercept=True)
regTreat = LinearRegression(fit_intercept=True)

# Fit regression line
regCont.fit( perc_filtered['mean'].values.reshape(-1, 1), perc_filtered['D1_Cont'].values.reshape(-1, 1))
regTreat.fit(perc_filtered['mean'].values.reshape(-1, 1), perc_filtered['D1_Treat'].values.reshape(-1, 1))
print(regTreat.coef_)
print(regCont.coef_)

# Correct Bias
# df['D1_Cont'] = df['D1_Cont'] - regCont.intercept_
# df['D1_Treat'] = df['D1_Treat'] - regTreat.intercept_

# Compute theta in radians
theta = calc_theta(regTreat.coef_[0][0], regCont.coef_[0][0])
print('theta', theta)
# theta=0.1

# Compute rotation Matrix
rotation_matrix_cont = compute_rot_mat(theta, regCont.coef_[0][0])
print(rotation_matrix_cont)

# rotate data
new_cont = rotate(rotation_matrix_cont, df[['mean', 'D1_Cont']])

# reassign data
new_data = pd.DataFrame({'Cont': new_cont[:, 1], 'Treat': df['D1_Treat']})
new_data.loc[:, 'mean'] = new_data.mean(axis=1)

fig, ax = plt.subplots()
new_data.plot('mean', 'Cont', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red', alpha=0.95, s=3)
new_data.plot('mean', 'Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='blue', alpha=0.25, s=15)

plt.plot([0, 5000], [0.0, regCont.coef_ * 5000], linestyle='--', color='white')
plt.plot([0, 5000], [0.0, regTreat.coef_ * 5000], linestyle='--', color='white');

ax.plot([0, 10000], [0, 10000], color='black')
ax.set_title("TMM Normalized, rotated");




d1.loc[:, 'diff'] = new_data['Cont'].sub(new_data['Treat']).abs()




d = df.copy()




d.head()




d.loc[:, 'diff'] = d['D1_Cont'].sub(d['D1_Treat']).abs()




d.describe(percentiles=np.linspace(0, 0.99, num=9))




d.describe(percentiles=np.linspace(0, 0.99, num=9)).drop(['count', 'mean', 'std', 'min', 'max']).plot(style=['--', '-', '--', '-']);




new_data[(new_data.apply(lambda x: x < 5000)) & (new_data.apply(lambda x: x > 100))].plot(kind='hist', bins=50, alpha=0.2)




euclidean_distance




new_df.dropna(axis=0, inplace=False)




# This function returns the TMM normalized of ALL sample. The next function needs to normalize by paired time point.
def TMM_all_samples(raw_df, 
                    threshold=10,
                    logfold_change_cutoff = 0.3,
                    absolute_intensity_cutoff = 0.05):
    
    frame = raw_df.copy()
    original_df = raw_df.copy()
    
    # Remove genes that don't have expression, and save them as a separate df to be added in later.
    df = frame.loc[(frame >= threshold).all(axis=1)]
    non_expressed_genes = frame.drop(df.index)
    
    ## Calculate factor quantiles
    ### PERFORM TRIMMING ###
    sample_column_names = df.columns.tolist()
    
    ref_col = sample_column_names[0]
    non_ref_cols = sample_column_names[1:]
    
    # 1. Calculate reference column statistics -- This could be more sophistocated (PCA to find middle most sample?
    # n_kref = column sum for reference column (first column)    
    reference_column_sum = df.loc[:, ref_col].sum()
    reference_column = df[ref_col]

    
    # 2. Calculate the gene-wise logfold changes
    genewise_logfold_change = df.apply(lambda col: np.log2((col/col.sum())/(reference_column/reference_column_sum)))
    
    
    # 3. Calculate the genewise Absolute Intensity
    genewise_abs_intensity = df.apply(lambda col: (np.log2((col/col.sum())*(reference_column/reference_column_sum)))/2.)

    
    # 4 Trim logfold and absolute     
    logfold_mask = genewise_logfold_change.apply(lambda col: col > genewise_logfold_change.quantile(logfold_change_cutoff, axis=1))
    abs_instensity_mask = genewise_abs_intensity.apply(lambda col: col > genewise_abs_intensity.quantile(absolute_intensity_cutoff, axis=1))
    
    import pdb;pdb.set_trace()

    # Remove rows with masked values
    df = df[logfold_mask]
    df = df[abs_instensity_mask]
    df_clean = df.dropna(axis=1)
    
    
    #######################################################
    """ From here we have a trimmed data frame of values that we can use to calculate normalization factors for."""
    # w-r_gk
    weight_r_gk = lambda col: ((col.sum() - col) / (col.sum() * col)) + ((reference_column_sum - reference_column) / (reference_column_sum * reference_column))
    w_r_gk = df_clean.apply(weight_r_gk)
    
    # m-r_gk
    mean_r_gk = lambda col: (np.log2(col / col.sum()) ) / (np.log2(reference_column / reference_column_sum))
    m_r_gk = df_clean.apply(mean_r_gk)
    
    normalization_factors = np.sum(w_r_gk * m_r_gk) / np.sum(w_r_gk)
    
    normalized_df = frame.multiply(normalization_factors, axis=1)
    
    return normalized_df




TMM_all_samples(raw_df)




raw_df.loc[(raw_df > 10).any(axis=1)]




ser =[0.1, 0.5, 0]




raw_df.mul(ser, axis=1)




calcNormFactors <- function(object, ...)
UseMethod("calcNormFactors")

calcNormFactors.DGEList <- function(object, method=c("TMM","RLE","upperquartile","none"), refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, ...)
#	Scale normalization of RNA-Seq data, for DGEList objects
#	Created 2 October 2014.  Last modified 27 August 2015.
{
	object$samples$norm.factors <- calcNormFactors(object=object$counts, lib.size=object$samples$lib.size, method=method, refColumn=refColumn, logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff, p=p)
	object
}

calcNormFactors.default <- function(object, lib.size=NULL, method=c("TMM","RLE","upperquartile","none"), refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75, ...)
#	Scale normalization of RNA-Seq data, for count matrices
#	Mark Robinson.  Edits by Gordon Smyth.
#	Created October 22 October 2009 by Mark Robinson.
#	Last modified 31 July 2015.
{
#	Check object
	x <- as.matrix(object)
	if(any(is.na(x))) stop("NA counts not permitted")

#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(x)
	if(any(is.na(lib.size))) stop("NA lib.sizes not permitted")

#	Check method
	method <- match.arg(method)

#	Remove all zero rows
	allzero <- .rowSums(x>0, nrow(x), ncol(x)) == 0
	if(any(allzero)) x <- x[!allzero,,drop=FALSE]

#	Degenerate cases
	if(nrow(x)==0 || ncol(x)==1) method="none"

#	Calculate factors
	f <- switch(method,
		TMM = {
			f75 <- .calcFactorQuantile(data=x, lib.size=lib.size, p=0.75)
			if( is.null(refColumn) ) refColumn <- which.min(abs(f75-mean(f75)))
			if(length(refColumn)==0 | refColumn < 1 | refColumn > ncol(x)) refColumn <- 1
			f <- rep(NA,ncol(x))
			for(i in 1:ncol(x))
				f[i] <- .calcFactorWeighted(obs=x[,i],ref=x[,refColumn], libsize.obs=lib.size[i], libsize.ref=lib.size[refColumn], logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
			f
		},
		RLE = .calcFactorRLE(x)/lib.size,
		upperquartile = .calcFactorQuantile(x,lib.size,p=p),
		none = rep(1,ncol(x))
	)

#	Factors should multiple to one
	f <- f/exp(mean(log(f)))

#	Output
	f
}


.calcFactorRLE <- function (data)
{
	gm <- exp(rowMeans(log(data)))
	apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorQuantile <- function (data, lib.size, p=0.75)
{
#	i <- apply(data<=0,1,all)
#	if(any(i)) data <- data[!i,,drop=FALSE]
	y <- t(t(data)/lib.size)
	f <- apply(y,2,function(x) quantile(x,p=p))
#	f/exp(mean(log(f)))
}

.calcFactorWeighted <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
#	TMM between two libraries
{
	obs <- as.numeric(obs)
	ref <- as.numeric(ref)

	if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
	if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref

	logR <- log2((obs/nO)/(ref/nR))			# log ratio of expression, accounting for library size
	absE <- (log2(obs/nO) + log2(ref/nR))/2	# absolute expression
	v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref	 # estimated asymptotic variance

#	remove infinite values, cutoff based on A
	fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

	logR <- logR[fin]
	absE <- absE[fin]
	v <- v[fin]

	if(max(abs(logR)) < 1e-6) return(1)

#	taken from the original mean() function
	n <- length(logR)
	loL <- floor(n * logratioTrim) + 1
	hiL <- n + 1 - loL
	loS <- floor(n * sumTrim) + 1
	hiS <- n + 1 - loS

#	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
#	a fix from leonardo ivan almonacid cardenas, since rank() can return
#	non-integer values when there are a lot of ties
	keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

	if(doWeighting)
		f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
	else
		f <- mean(logR[keep], na.rm=TRUE)

#	Results will be missing if the two libraries share no features with positive counts
#	In this case, return unity
	if(is.na(f)) f <- 0
	2^f
}


