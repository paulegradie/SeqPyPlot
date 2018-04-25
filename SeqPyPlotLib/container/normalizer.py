"""Somone's implementation of TMM. Can take this directly from smyths code"""

import numpy as np


def compute_scaling_factor(obs, ref,
                           log_ratio_trim=0.3,
                           sum_trim=0.05,
                           weighting=True,
                           a_cutoff=-1e10):
    """
    This is largely taken from the Smyth implementation in R.
    """
    if all(obs == ref):
        return 1

    obs_sum = obs.sum()
    ref_sum = ref.sum()

    # log ration of expression accounting for library size
    lr = np.log2((obs / obs_sum) / (ref / ref_sum))
    
    # absolute expression
    ae = (np.log2(obs / obs_sum) + np.log2(ref / ref_sum)) / 2
    
    # estimated asymptotic variance
    v = (obs_sum - obs) / obs_sum / obs + (ref_sum - ref) / ref_sum / ref
    
    # create mask
    m = np.isfinite(lr) & np.isfinite(ae) & (ae > a_cutoff)
    
    # drop the masked values
    lr = lr[m]
    ae = ae[m]
    v = v[m]
    assert len(lr) == len(ae) == len(v)

    n = len(lr)
    lo_l = np.floor(n * log_ratio_trim) + 1
    hi_l = n + 1 - lo_l
    lo_s = np.floor(n * sum_trim) + 1
    hi_s = n + 1 - lo_s
    k = ((lr.rank(method="first") >= lo_l) & (lr.rank(method="first") <= hi_l)) \
          & ((ae.rank(method="first") >= lo_s) & (ae.rank(method="first") <= hi_s))

    if weighting:
        return 2**(sum(lr[k] / v[k]) / sum(1 / v[k]))
    else:
        return 2**(lr[k].mean())


def extract_usable_data(df):

    d = df.copy()
    d.fillna(0, inplace=True)
    
    nonzero_df = df.where(df > 0).dropna(axis=0)
    zero_df = df.where(df == 0).dropna(axis=0)

    return nonzero_df, zero_df


def norm_tmm(original_df,
             log_ratio_trim=0.3,
             absolute_intensity_trim=0.05,
             weighting=True,
             a_cutoff=-1e10):
    """
    Trimmed Mean of M-values (TMM) is the weighted mean of log ratios between
    this test and the reference, after exclusion of the most expressed genes
    and the genes with the largest log ratios.

    Parameters:
    -----------
    - original_df: experiment object.
    - ref_col: reference column from which to scale others.
    - log_ratio_trim: amount of trim to use on log-ratios.
    - absolute_intensity_trim : amount of trim to use on combined absolute values.
    - weighting: whether to compute weights.
    - a_cutoff: cutoff of absolute expression values.
    """
    
    df = original_df.copy()
    nonzero_df, zero_df = extract_usable_data(df)

    # Always use control as ref
    ref_col = df.columns[0]

    kwargs = {"ref": nonzero_df[ref_col],
              "log_ratio_trim": log_ratio_trim,
              "sum_trim": absolute_intensity_trim,
              "weighting": weighting,
              "a_cutoff": a_cutoff}

    sf_tmm = nonzero_df.apply(compute_scaling_factor, **kwargs)
    
    remerged_df = nonzero_df.merge(zero_df, how='right')

    # apply scaling
    df = df.div(sf_tmm, axis=1)
    return df














# def _scaling_factor(df, quartile=0.60):
#     """
#     Parameters:
#     -----------
#     - df: zeroed rows removed
#     - quartile: percent used to claculate trimming of outer quartile data
#     """
#     library_sizes = df.sum()

#     result_df = df.div(library_sizes, axis=1)
#     result_df = result_df.quantile(quartile)

#     # factors multiple to one
#     scaling_factor = result_df.div(np.exp(np.mean(np.log(result_df))))
#     return scaling_factor





# def _scaling_factor_for_tmm(df, quartile=0.75):
#     """
#     Parameters:
#     -----------
#     - df: zeroed rows removed
#     - q: quartile
#     """
#     lib_size = df.sum()
#     y = df.div(lib_size, axis=0)
#     # fill nans with 0
#     y = y.dropna(how="all").fillna(0)
#     y = y.quantile(quartile)
#     # factors multiple to one
#     scaling_factor = y.div(np.exp(np.mean(np.log(y))))
#     return scaling_factor

# def norm_tmm(exp_obj,
#              ref_col=None,
#              log_ratio_trim=0.3,
#              sum_trim=0.05,
#              weighting=True,
#              a_cutoff=-1e10):
#     """
#     Trimmed Mean of M-values (TMM) is the weighted mean of log ratios between
#     this test and the reference, after exclusion of the most expressed genes
#     and the genes with the largest log ratios.
    
#     Parameters:
#     -----------
#     - exp_obj: experiment object.
#     - ref_col: reference column from which to scale others.
#     - log_ratio_trim: amount of trim to use on log-ratios.
#     - sum_trim: amount of trim to use on combined absolute values.
#     - weighting: whether to compute weights.
#     - a_cutoff: cutoff of absolute expression values.
#     """
#     df = exp_obj.counts_df.copy()
#     # remove zeros
#     nz = df.where(df > 0)
#     nz = nz.dropna(how="all").fillna(0)
#     # reference column
#     if ref_col is None:
#         # quantile factors
#         sf_q = _sf_q(nz)
#         ref_col = (abs(sf_q - np.mean(sf_q))).idxmin()
#     # try:
#     kwargs = {"ref": nz[ref_col],
#               "log_ratio_trim": log_ratio_trim,
#               "sum_trim": sum_trim,
#               "weighting": weighting,
#               "a_cutoff": a_cutoff}
#     # except KeyError:
#         # revert back to auto?
#     sf_tmm = nz.apply(_sf_tmm, **kwargs)
#     # apply scaling
#     df = df.div(sf_tmm, axis=1)
#     return df

# def norm_q(exp_obj, q=0.75):
#     """
#     Ported from edgeR and still needs to be validated. Also, maybe compare
#     edgeR method to limma implementation.
    
#     Quantile normalization.
    
#     Parameters:
#     -----------
#     - exp_obj: experiment object.
#     - q: quantile.
#     """
#     df = exp_obj.counts_df.copy()
#     # remove zeros
#     nz = df.where(df > 0)
#     nz = nz.dropna(how="all").fillna(0)
#     sf_q = _sf_q(nz, q)
#     # apply scaling
#     df = df.div(sf_q, axis=1)
#     return df