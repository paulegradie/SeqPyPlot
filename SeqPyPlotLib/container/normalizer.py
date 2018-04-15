"""Somone's implementation of TMM. Can take this directly from smyths code"""

import numpy as np


def _scaling_factor_for_tmm(obs, ref,
                            log_ratio_trim=0.3,
                            sum_trim=0.05,
                            weighting=True,
                            a_cutoff=-1e10):
    """
    Called by `norm_tmm`.
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


def _scaling_factor(df, q=0.75):
    """
    Parameters:
    -----------
    - df: zeroed rows removed
    - q: quartile
    """
    # Claculate column sums
    library_sizes = df.sum()

    # Divide rows by column sums and return to
    # original orientation
    result = df.T.div(library_sizes, axis=0).T

    # # fill nans with 0
    # result = result.dropna(how="all").fillna(0)
    result = result.quantile(q)

    # factors multiple to one
    scaling_factor = result.div(np.exp(np.mean(np.log(result))))
    return scaling_factor


def extract_usable_data(df):

    nonzero_df = df[df.values.sum(axis=1) != 0]
    zero_df = df[df.values.sum(axis=1) == 0]

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

    #Take whichever column has the smallest
    ref_col = df.columns[0]

    # try:
    kwargs = {"ref": nonzero_df[ref_col],
              "log_ratio_trim": log_ratio_trim,
              "sum_trim": absolute_intensity_trim,
              "weighting": weighting,
              "a_cutoff": a_cutoff}

    sf_tmm = nonzero_df.apply(_scaling_factor_for_tmm, **kwargs)
    # apply scaling
    df = df.div(sf_tmm, axis=1)
    return df
