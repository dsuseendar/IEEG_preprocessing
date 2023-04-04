from functools import singledispatch

import numpy as np
from mne.utils import logger, verbose
from mne.epochs import BaseEpochs
from mne import Epochs

from tqdm import tqdm
from numba import njit
from ..utils.utils import parallelize


def sum_squared(x: np.ndarray) -> np.ndarray:
    """Compute norm of an array.

    Parameters
    ----------
    x : array
        Data whose norm must be found.
    Returns
    -------
    value : float
        Sum of squares of the input array x.
    """
    x_flat = x.ravel()
    return np.dot(x_flat, x_flat)


# @njit([t.Tuple((float64[:, ::1], complex128[:, ::1]))(
#     float64[:, :], complex128[:, :, ::1])], nogil=True, boundscheck=True,
#     fastmath=True, cache=True)
def sine_f_test(window_fun: np.ndarray, x_p: np.ndarray
                ) -> (np.ndarray, np.ndarray):
    """computes the F-statistic for sine wave in locally-white noise.

    Parameters
    ----------
    window_fun : array
        The tapers used to calculate the multitaper spectrum.
    x_p : array
        The tapered time series.

    Returns
    -------
    f_stat : array
        The F-statistic for each frequency.
    A : array
        The amplitude of the sine wave at each frequency.
    """
    # drop the even tapers
    n_tapers = len(window_fun)
    tapers_odd = np.arange(0, n_tapers, 2)
    tapers_even = np.arange(1, n_tapers, 2)
    tapers_use = window_fun[tapers_odd]

    # sum tapers for (used) odd prolates across time (n_tapers, 1)
    H0 = np.sum(tapers_use, axis=1)

    # sum of squares across tapers (1, )
    H0_sq = sum_squared(H0)

    # sum of the product of x_p and H0 across tapers (1, n_freqs)
    exp_H0 = np.reshape(H0, (1, -1, 1))
    x_p_H0 = np.sum(x_p[:, tapers_odd, :] * exp_H0, axis=1)

    # resulting calculated amplitudes for all freqs
    A = x_p_H0 / H0_sq

    # figure out which freqs to remove using F stat

    # estimated coefficient
    x_hat = A * np.reshape(H0, (-1, 1))

    # numerator for F-statistic
    num = (n_tapers - 1) * (A * A.conj()).real * H0_sq
    # denominator for F-statistic
    den = (np.sum(np.abs(x_p[:, tapers_odd, :] - x_hat) ** 2, 1) +
           np.sum(np.abs(x_p[:, tapers_even, :]) ** 2, 1))
    den = np.where(den == 0, np.inf, den)
    f_stat = num / den

    return f_stat, A


def _log_rescale(baseline, mode='mean'):
    """Log the rescaling method."""
    if baseline is not None:
        msg = 'Applying baseline correction (mode: %s)' % mode
    else:
        msg = 'No baseline correction applied'
    return msg


@singledispatch
def rescale(data: np.ndarray, basedata: np.ndarray, mode: str = 'mean',
            copy: bool = False) -> np.ndarray:
    """Rescale (baseline correct) data.

    Parameters
    ----------
    data : array
        It can be of any shape. The only constraint is that the last
        dimension should be time.
    basedata : array
        It can be of any shape. The last dimension should be time, and the
        other dimensions should be the same as data.
    mode : 'mean' | 'ratio' | 'logratio' | 'percent' | 'zscore' | 'zlogratio',\
        default 'mean', optional
        Perform baseline correction by
        - subtracting the mean of baseline values ('mean')
        - dividing by the mean of baseline values ('ratio')
        - dividing by the mean of baseline values and taking the log
        ('logratio')
        - subtracting the mean of baseline values followed by dividing by
        the mean of baseline values ('percent')
        - subtracting the mean of baseline values and dividing by the
        standard deviation of baseline values ('zscore')
        - dividing by the mean of baseline values, taking the log, and
        dividing by the standard deviation of log baseline values ('zlogratio')
    copy : bool, optional
        Whether to return a new instance or modify in place.

    Returns
    -------
    data_scaled: array
        Array of same shape as data after rescaling.
    """
    if copy:
        data = data.copy()

    match mode:
        case 'mean':
            def fun(d, m, s):
                d -= m
        case 'ratio':
            def fun(d, m, s):
                d /= m
        case 'logratio':
            def fun(d, m, s):
                d /= m
                np.log10(d, out=d)
        case 'percent':
            def fun(d, m, s):
                d -= m
                d /= m
        case 'zscore':
            def fun(d, m, s):
                d -= m
                d /= s
        case 'zlogratio':
            def fun(d, m, s):
                d /= m
                np.log10(d, out=d)
                d /= s
        case _:
            raise NotImplementedError()
    mean = np.mean(basedata, axis=-1, keepdims=True)
    std = np.std(basedata, axis=-1, keepdims=True)
    fun(data, mean, std)
    return data


@rescale.register
@verbose
def _(line: BaseEpochs, baseline: BaseEpochs, mode: str = 'mean',
      copy: bool = False, picks: list = 'data', verbose=None) -> Epochs:
    """Rescale (baseline correct) Epochs"""
    if copy:
        line: Epochs = line.copy()
    if verbose is not False:
        msg = _log_rescale(baseline, mode)
        logger.info(msg)

    # Average the baseline across epochs
    basedata = np.mean(baseline.pick(picks)._data, axis=0, keepdims=True)
    line.pick(picks)._data = rescale(line.pick(picks)._data, basedata, mode,
                                     False)
    return line


def make_data_same(data_fix: np.ndarray, data_like: np.ndarray,
                   ignore_axis: int | tuple[int] = None) -> np.ndarray:
    """Make data_fix the same shape as data_like.

    Reshapes data_fix to match the shape of data_like, ignoring the
    specified axis. If data_fix is already the same shape as data_like,
    it is returned unchanged. If data_fix is larger in the ignored
    axis and smaller in another unignored axis, the data is split along
    the ignored axis and appended to the end of the unignored one.

    Parameters
    ----------
    data_fix : array
        The data to reshape.
    data_like : array
        The data to match the shape of.
    ignore_axis : int | tuple[int], optional
        The axis to ignore. If None, all axes are used.

    Returns
    -------
    data_fix : array
        The reshaped data.
    """
    if ignore_axis is None:
        ignore_axis = tuple()
    elif isinstance(ignore_axis, int):
        ignore_axis = (ignore_axis,)
    else:
        ignore_axis = tuple(ignore_axis)

    data = data_fix.copy()
    shape_fix = list(data_fix.shape)
    shape_like = data_like.shape
    if len(shape_fix) != len(shape_like):
        raise ValueError('data_fix and data_like must have the same number of '
                         'dimensions')
    for i, (s1, s2) in enumerate(zip(shape_fix, shape_like)):
        if i not in ignore_axis and s1 != s2:
            if len(ignore_axis) == 0:
                raise ValueError('data_fix and data_like must have the same '
                                 'shape if ignore_axis is None')
            repeats = int(s2 / s1)
            reduce_axis = ignore_axis[0]
            while shape_fix[reduce_axis] % repeats != 0:
                repeats -= 1
            shape_fix[reduce_axis] = shape_fix[reduce_axis] // repeats
            shape_fix[i] = shape_fix[i] * repeats
            data = np.reshape(data, shape_fix)
            out_pad = [(0, 0) if j != i else (0, s2 - (s1 * repeats)) for j in
                       range(len(shape_fix))]
            data = np.pad(data, out_pad, 'wrap')

    return data


@njit(fastmath=True)
def mean(data: np.ndarray, axis: int = 0) -> np.ndarray:
    """Calculate the mean of an array using numba compatable methods.

    Parameters
    ----------
    data : array
        The data to calculate the mean of.
    axis : int, optional
        The axis to calculate the mean across.

    Returns
    -------
    mean : array
        The mean of the data.
    """
    return np.sum(data, axis=axis) / data.shape[axis]


def time_perm_cluster(sig1: np.ndarray, sig2: np.ndarray, n_perm: int = 1000,
                      tails: int = 1, axis: int = 0) -> np.ndarray:
    
    if len(sig1.shape) > 2:
        if any(sig1.shape[1:-2] != sig2.shape[1:-2]):
            trials = int(sig2[:, 0, ..., :].size / sig1.shape[-1])
            sigB = np.full((trials, *sig1.shape[1:]), np.nan)
            for i in range(sig1.shape[1]):
                sigB[:, i, :].flat = sig2[:, i, :].flat
            sig2 = sigB.copy()
            del trials, sigB


def time_perm_shuffle(sig1: np.ndarray, sig2: np.ndarray, n_perm: int = 1000,
                      tails: int = 1, axis: int = 0) -> np.ndarray:
    """Time permutation cluster test between two time series.

    The test is performed by shuffling the trials of the two time series and
    calculating the difference between the two groups at each time point. The
    p-value is the proportion of times the difference between the two groups
    is greater than the original observed difference. The number of trials in
    each group does not need to be the same.

    Parameters
    ----------
    sig1 : array, shape (trials, ..., time)
        Active signal. The first dimension is assumed to be the trials
    sig2 : array, shape (trials, ..., time)
        Passive signal. The first dimension is assumed to be the trials
    n_perm : int, optional
        The number of permutations to perform.
    tails : int, optional
        The number of tails to use. 1 for one-tailed, 2 for two-tailed.

    Returns
    -------
    p : np.ndarray, shape (..., time)
        The p-values for each time point.
        """
    # Concatenate the two signals for trial shuffling
    all_trial = np.concatenate((sig1, sig2), axis=axis)
    labels = np.concatenate((np.full(sig1.shape[axis], False, dtype=bool),
                             np.full(sig2.shape[axis], True, dtype=bool)))

    perm_labels = np.array(
        [np.random.permutation(labels) for _ in range(n_perm)])

    # Calculate the observed difference
    obs_diff = np.mean(sig1, axis) - np.mean(sig2, axis)

    # Calculate the difference between the two groups averaged across
    # trials at each time point
    fake_sig1 = np.array([np.take(all_trial, np.where(  # where False
        np.invert(perm_labels[i]))[0], axis=axis) for i in range(n_perm)])
    fake_sig2 = np.array([np.take(all_trial, np.where(  # where True
        perm_labels[i])[0], axis=axis) for i in range(n_perm)])
    diff = np.mean(fake_sig1, axis=axis + 1) - np.mean(fake_sig2,
                                                       axis=axis + 1)
    if tails == 1:
        larger = diff > obs_diff
    elif tails == 2:
        larger = np.abs(diff) > np.abs(obs_diff)
    else:
        raise ValueError('tails must be 1 or 2')

    # Calculate the p-value
    p = np.mean(larger, axis=0)

    return p


def time_perm_shuff2(sig1: np.ndarray, sig2: np.ndarray, n_perm: int = 1000,
                      tails: int = 1, axis: int = 0, verbose: bool = True) -> np.ndarray:
    # Concatenate the two signals for trial shuffling
    all_trial = np.concatenate((sig1, sig2), axis=axis)
    labels = np.concatenate((np.zeros(sig1.shape[axis]),
                             np.ones(sig2.shape[axis])))

    all_labels = [np.random.permutation(labels) for _ in range(n_perm)]
    if verbose:
        all_labels = tqdm(all_labels)

    # Calculate the observed difference
    obs_diff = np.mean(sig1, axis) - np.mean(sig2, axis)
    # larger = np.full((n_perm,) + tuple(obs_diff.shape), False)

    larger = parallelize(sample_test, all_labels, all_trial=all_trial,
                         obs_diff=obs_diff, tails=tails, axis=axis, n_jobs=-1)

    return np.mean(np.array(larger), axis=0)


def sample_test(labels: np.ndarray, all_trial: np.ndarray, obs_diff: np.ndarray,
                 tails: int = 1, axis: int = 0):
    # Calculate the difference between the two groups averaged across
    # trials at each time point
    fake_sig1 = np.take(all_trial, np.where(labels == 0)[axis], axis=axis)
    fake_sig2 = np.take(all_trial, np.where(labels == 1)[axis], axis=axis)
    diff = mean(fake_sig1, axis=axis) - mean(fake_sig2, axis=axis)
    if tails == 1:
        return diff > obs_diff
    elif tails == 2:
        return np.abs(diff) > np.abs(obs_diff)
    else:
        raise ValueError('tails must be 1 or 2')
