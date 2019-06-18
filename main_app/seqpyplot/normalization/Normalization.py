
# coding: utf-8

# # Heteroskedasticity correction notebook using linear transformations rotation matrices
#
# This notebook works through an idea for correcting heteroskedasticity which uses rotation matrices to align samples.
#
# For a given sample pair, the mean across all samples is computed. Then for each sample and the mean, linear regression is used to compute a line of best fit through the sample. The result is a coefficient for each sample. If these coeffieicients do not match, we can use a rotation matrix to rotate the data about the origin until they match, resulting in proper alignment between the control and treated samples.
#
# Further correction can be introduced by first negating any bias computed during linear regression.



import pandas as pd
import numpy as np
from random import randint as rand
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

from sklearn.linear_model import LinearRegression

from scipy.linalg import svd

from seqpyplot.container.data_container import DataContainer
from seqpyplot.parsers.config_parser import config_parser
from pathlib import Path

from matplotlib import rcParams

rcParams['figure.figsize'] = (15, 15)


# ## Example
# Example of rotating a line around the origin to match the slope of another line



slope1 = 1.1
slope2 = 2.0

line1 = np.array([slope1 * x for x in range(10)])
line2 = np.array([slope2 * x for x in range(10)])

xs = list(range(10))




plt.plot(xs, line1, color='black');
plt.plot(xs, line2, color='red');




def calc_theta(coef1, coef2):
    "returns radians"
    return np.abs(
       np.arctan(np.abs(coef1 - coef2) / (1. + (coef1 * coef2)))
   )




angle_diff = calc_theta(slope1, slope2)
angle_diff




def compute_rot_mat(rad, coef=.5):
    if coef < 1.0:
        rotation_matrix = np.array([[np.cos(rad), -np.sin(rad)],
                                    [np.sin(rad), np.cos(rad)]])
    else:
        rotation_matrix = np.array([[np.cos(rad), np.sin(rad)],
                                    [-np.sin(rad), np.cos(rad)]])
    return rotation_matrix




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


# #### Demonstration of aligning lines using results from rotation matrix



plt.scatter(xs, line1, color='red')
plt.scatter(new_line1[:, 0], new_line1[:, 1], color='green')

plt.scatter(xs, line2, color='black', alpha=0.4, s=95);
for (x1, y1), (x2, y2) in zip(zip(xs, line1), new_line1):

    plt.plot([x1, x2], [y1, y2], linestyle='--', color='y')


# # Examle using actual Data No TMM
#
# 1. Correct bias
# 2. Rotate data



config = 'examples/example_config.ini'
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
d1.plot('mean', 'D1_Cont', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='blue')
d2.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red');
ax.set_title("Raw unnormalized data");




regCont = LinearRegression(fit_intercept=True)
regTreat = LinearRegression(fit_intercept=True)




regCont.fit(d1['D1_Cont'].values.reshape(-1, 1), d1['mean'].values.reshape(-1, 1))
regTreat.fit(d2['D1_Treat'].values.reshape(-1, 1), d1['mean'].values.reshape(-1, 1))




regCont.coef_, regCont.intercept_




regTreat.coef_, regTreat.intercept_


# #### Correct Bias



d1['D1_Cont'] = d1['D1_Cont'] - regCont.intercept_
d2['D1_Treat'] = d2['D1_Treat'] - regTreat.intercept_




fig, ax = plt.subplots()
d1.plot('mean', 'D1_Cont', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='blue')
d2.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red')
plt.plot([0, 5000], [0.0, regCont.coef_ * 5000], linestyle='--', color='black')
plt.plot([0, 5000], [0.0, regTreat.coef_ * 5000], linestyle='--', color='black');
ax.set_title("bias corrected, with best fit lines");


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


# ### Realigned Data



fig, ax = plt.subplots()
d1.plot('mean', 'D1_Cont', kind='scatter', xlim=(0, 20000), ylim=(0, 20000), ax=ax, color='blue')
d2_cor.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red');
ax.set_title("No TMM, linearly transformed");


# # Use TMM before linear transformation



# load the data container_obj
container_obj = DataContainer(config_obj)
data, ercc_data = container_obj.parse_input()

data = container_obj.normalize_file_pairs(data) # Single df of normalized data

data.head()

df = data[['D1_Cont', 'D1_Treat']]
df.loc[:, 'mean'] = df.mean(axis=1)

d1 = df[['D1_Cont', 'mean']]
d2 = df[['D1_Treat', 'mean']]




fig, ax = plt.subplots()
d1.plot('mean', 'D1_Cont', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='blue')
d2.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red');
ax.set_title("TMM Normalized");




regCont = LinearRegression(fit_intercept=True)
regTreat = LinearRegression(fit_intercept=True)

regCont.fit(d1['D1_Cont'].values.reshape(-1, 1), d1['mean'].values.reshape(-1, 1))
regTreat.fit(d2['D1_Treat'].values.reshape(-1, 1), d1['mean'].values.reshape(-1, 1))




regCont.coef_, regCont.intercept_




regTreat.coef_, regTreat.intercept_




d2.head()




d1['D1_Cont'] = d1['D1_Cont'] - regCont.intercept_
d2['D1_Treat'] = d2['D1_Treat'] - regTreat.intercept_




fig, ax = plt.subplots()
d1.plot('mean', 'D1_Cont', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='blue')
d2.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red')
plt.plot([0, 5000], [0.0, regCont.coef_ * 5000], linestyle='--', color='black')
plt.plot([0, 5000], [0.0, regTreat.coef_ * 5000], linestyle='--', color='white');
ax.set_title("bias corrected, with best fit lines");


# Compute linear transformation



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
d1.plot('mean', 'D1_Cont', kind='scatter', xlim=(0, 20000), ylim=(0, 20000), ax=ax, color='blue')
d2_cor.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red');
ax.set_title("With TMM, linearly transformed");




d2_cor.head()




reg = LinearRegression(fit_intercept=True)
reg.fit(d2_cor['mean'].values.reshape(-1, 1), d2_cor['D1_Treat'].values.reshape(-1, 1), )




fig, ax = plt.subplots()
d1.plot('mean', 'D1_Cont', kind='scatter', xlim=(0, 20000), ylim=(0, 20000), ax=ax, color='blue')
d2_cor.plot('mean', 'D1_Treat', kind='scatter', xlim=(0, 5000), ylim=(0, 5000), ax=ax, color='red');
plt.plot([0, 5000], [0.0, reg.coef_ * 5000], linestyle='--', color='white')
# plt.plot([0, 5000], [0.0, regCont.coef_ * 5000], linestyle='--', color='black');
# plt.plot([0, 5000], [0.0, regTreat.coef_ * 5000], linestyle='--', color='blue');

ax.set_title("With TMM, linearly transformed");




calc_theta(1.0, regCont.coef_)




calc_theta(1.0, regTreat.coef_)


# The most straightforward approach is probably going to be to transform the data so that each one line is brought to the slope of another. Once both samples are transformed, then the means will be dropped, and then recomputed. This will cause the data to be reordered (since we order based on the mean). Once plotted, the redistributed data should then be in a form that is ready for filtering.

# # Choosing thresholds for outlier cutoffs
#
# Before performing rotation correction, its a good idea to eliminate some of the worst outlier offenders before performing linear regression. To do this, we can calculate summary statistics on the sample pair, and



config = 'examples/example_config.ini'
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
df['D1_Cont'] = df['D1_Cont'] - regCont.intercept_
df['D1_Treat'] = df['D1_Treat'] - regTreat.intercept_

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


