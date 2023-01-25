import operator
import gc
import os.path as op
import pandas as pd
from os import PathLike as PL
from os import environ
from typing import List, TypeVar, Iterable, Union

import matplotlib as mpl
from matplotlib.pyplot import Axes
from mne.io import Raw
from mne.utils import config, logger, verbose
import numpy as np
from joblib import Parallel, delayed, cpu_count
from tqdm import tqdm

try:
    mpl.use("TkAgg")
except ImportError:
    pass

HOME = op.expanduser("~")
LAB_root = op.join(HOME, "Box", "CoganLab")
PathLike = TypeVar("PathLike", str, PL)


# plotting funcs

def figure_compare(raw: List[Raw], labels: List[str], avg: bool = True,
                   n_jobs: int = None, **kwargs):
    """Plots the psd of a list of raw objects"""
    if n_jobs is None:
        n_jobs = cpu_count() - 2
    for title, data in zip(labels, raw):
        title: str
        data: Raw
        psd = data.compute_psd(n_jobs=n_jobs, **kwargs,
                               n_fft=int(data.info['sfreq']))
        fig = psd.plot(average=avg, spatial_colors=avg)
        fig.subplots_adjust(top=0.85)
        fig.suptitle('{}filtered'.format(title), size='xx-large',
                     weight='bold')
        add_arrows(fig.axes[:2])
        gc.collect()


def add_arrows(axes: Axes):
    """add some arrows at 60 Hz and its harmonics"""
    for ax in axes:
        freqs = ax.lines[-1].get_xdata()
        psds = ax.lines[-1].get_ydata()
        for freq in (60, 120, 180, 240):
            idx = np.searchsorted(freqs, freq)
            # get ymax of a small region around the freq. of interest
            y = psds[(idx - 4):(idx + 5)].max()
            ax.arrow(x=freqs[idx], y=y + 18, dx=0, dy=-12, color='red',
                     width=0.1, head_width=3, length_includes_head=True)


def ensure_int(x, name='unknown', must_be='an int', *, extra=''):
    """Ensure a variable is an integer."""
    # This is preferred over numbers.Integral, see:
    # https://github.com/scipy/scipy/pull/7351#issuecomment-299713159
    extra = f' {extra}' if extra else extra
    try:
        # someone passing True/False is much more likely to be an error than
        # intentional usage
        if isinstance(x, bool):
            raise TypeError()
        x = int(operator.index(x))
    except TypeError:
        raise TypeError(f'{name} must be {must_be}{extra}, got {type(x)}')
    return x


def validate_type(item, types):
    try:
        if isinstance(types, TypeVar):
            check = isinstance(item, types.__constraints__)
        elif types is int:
            ensure_int(item)
            check = True
        elif types is float:
            check = is_number(item)
        else:
            check = isinstance(item, types)
    except TypeError:
        check = False
    if not check:
        raise TypeError(
            f"must be an instance of {types}, "
            f"got {type(item)} instead.")


def is_number(s) -> bool:
    if isinstance(s, str):
        try:
            float(s)
            return True
        except ValueError:
            return False
    elif isinstance(s, (np.number, int, float)):
        return True
    elif isinstance(s, pd.DataFrame):
        try:
            s.astype(float)
            return True
        except Exception:
            return False
    elif isinstance(s, pd.Series):
        try:
            pd.to_numeric(s)
            return True
        except Exception:
            return False
    else:
        return False


def parallelize(func: object, par_var: Iterable, n_jobs: int = None, *args,
                **kwargs) -> list:
    if n_jobs is None:
        n_jobs = cpu_count()
    elif n_jobs == -1:
        n_jobs = cpu_count()
    settings = dict(verbose=5,  # prefer='threads',
                    pre_dispatch=n_jobs)
    env = dict(**environ)
    if config.get_config('MNE_CACHE_DIR') is not None:
        settings['temp_folder'] = config.get_config('MNE_CACHE_DIR')
    elif 'TEMP' in env.keys():
        settings['temp_folder'] = env['TEMP']

    if config.get_config('MNE_MEMMAP_MIN_SIZE') is not None:
        settings['max_nbytes'] = config.get_config('MNE_MEMMAP_MIN_SIZE')
    else:
        settings['max_nbytes'] = get_mem()

    data_new = Parallel(n_jobs, **settings)(delayed(func)(
        x_, *args, **kwargs)for x_ in tqdm(par_var))
    return data_new


def get_mem() -> Union[float, int]:
    from psutil import virtual_memory
    ram_per = virtual_memory()[3]/cpu_count()
    return ram_per


###############################################################################
# Constant overlap-add processing class


def _check_store(store):
    if isinstance(store, np.ndarray):
        store = [store]
    if isinstance(store, (list, tuple)) and all(isinstance(s, np.ndarray)
                                                for s in store):
        store = _Storer(*store)
    if not callable(store):
        raise TypeError('store must be callable, got type %s'
                        % (type(store),))
    return store


class _COLA:
    r"""Constant overlap-add processing helper.
    Parameters
    ----------
    process : callable
        A function that takes a chunk of input data with shape
        ``(n_channels, n_samples)`` and processes it.
    store : callable | ndarray
        A function that takes a completed chunk of output data.
        Can also be an ``ndarray``, in which case it is treated as the
        output data in which to store the results.
    n_total : int
        The total number of samples.
    n_samples : int
        The number of samples per window.
    n_overlap : int
        The overlap between windows.
    window : str
        The window to use. Default is "hann".
    tol : float
        The tolerance for COLA checking.
    Notes
    -----
    This will process data using overlapping windows to achieve a constant
    output value. For example, for ``n_total=27``, ``n_samples=10``,
    ``n_overlap=5`` and ``window='triang'``::
        1 _____               _______
          |    \   /\   /\   /
          |     \ /  \ /  \ /
          |      x    x    x
          |     / \  / \  / \
          |    /   \/   \/   \
        0 +----|----|----|----|----|-
          0    5   10   15   20   25
    This produces four windows: the first three are the requested length
    (10 samples) and the last one is longer (12 samples). The first and last
    window are asymmetric.
    """

    @verbose
    def __init__(self, process, store, n_total, n_samples, n_overlap,
                 sfreq, window='hann', tol=1e-10, *, verbose=None):
        from scipy.signal import get_window
        n_samples = _ensure_int(n_samples, 'n_samples')
        n_overlap = _ensure_int(n_overlap, 'n_overlap')
        n_total = _ensure_int(n_total, 'n_total')
        if n_samples <= 0:
            raise ValueError('n_samples must be > 0, got %s' % (n_samples,))
        if n_overlap < 0:
            raise ValueError('n_overlap must be >= 0, got %s' % (n_overlap,))
        if n_total < 0:
            raise ValueError('n_total must be >= 0, got %s' % (n_total,))
        self._n_samples = int(n_samples)
        self._n_overlap = int(n_overlap)
        del n_samples, n_overlap
        if n_total < self._n_samples:
            raise ValueError('Number of samples per window (%d) must be at '
                             'most the total number of samples (%s)'
                             % (self._n_samples, n_total))
        if not callable(process):
            raise TypeError('process must be callable, got type %s'
                            % (type(process),))
        self._process = process
        self._step = self._n_samples - self._n_overlap
        self._store = _check_store(store)
        self._idx = 0
        self._in_buffers = self._out_buffers = None

        # Create our window boundaries
        window_name = window if isinstance(window, str) else 'custom'
        self._window = get_window(window, self._n_samples,
                                  fftbins=(self._n_samples - 1) % 2)
        self._window /= _check_cola(self._window, self._n_samples, self._step,
                                    window_name, tol=tol)
        self.starts = np.arange(0, n_total - self._n_samples + 1, self._step)
        self.stops = self.starts + self._n_samples
        delta = n_total - self.stops[-1]
        self.stops[-1] = n_total
        sfreq = float(sfreq)
        pl = 's' if len(self.starts) != 1 else ''
        logger.info('    Processing %4d data chunk%s of (at least) %0.1f sec '
                    'with %0.1f sec overlap and %s windowing'
                    % (len(self.starts), pl, self._n_samples / sfreq,
                       self._n_overlap / sfreq, window_name))
        del window, window_name
        if delta > 0:
            logger.info('    The final %0.3f sec will be lumped into the '
                        'final window' % (delta / sfreq,))

    @property
    def _in_offset(self):
        """Compute from current processing window start and buffer len."""
        return self.starts[self._idx] + self._in_buffers[0].shape[-1]

    @verbose
    def feed(self, *datas, verbose=None, **kwargs):
        """Pass in a chunk of data."""
        # Append to our input buffer
        if self._in_buffers is None:
            self._in_buffers = [None] * len(datas)
        if len(datas) != len(self._in_buffers):
            raise ValueError('Got %d array(s), needed %d'
                             % (len(datas), len(self._in_buffers)))
        for di, data in enumerate(datas):
            if not isinstance(data, np.ndarray) or data.ndim < 1:
                raise TypeError('data entry %d must be an 2D ndarray, got %s'
                                % (di, type(data),))
            if self._in_buffers[di] is None:
                # In practice, users can give large chunks, so we use
                # dynamic allocation of the in buffer. We could save some
                # memory allocation by only ever processing max_len at once,
                # but this would increase code complexity.
                self._in_buffers[di] = np.empty(
                    data.shape[:-1] + (0,), data.dtype)
            if data.shape[:-1] != self._in_buffers[di].shape[:-1] or \
                    self._in_buffers[di].dtype != data.dtype:
                raise TypeError('data must dtype %s and shape[:-1]==%s, '
                                'got dtype %s shape[:-1]=%s'
                                % (self._in_buffers[di].dtype,
                                   self._in_buffers[di].shape[:-1],
                                   data.dtype, data.shape[:-1]))
            logger.debug('    + Appending %d->%d'
                         % (self._in_offset, self._in_offset + data.shape[-1]))
            self._in_buffers[di] = np.concatenate(
                [self._in_buffers[di], data], -1)
            if self._in_offset > self.stops[-1]:
                raise ValueError('data (shape %s) exceeded expected total '
                                 'buffer size (%s > %s)'
                                 % (data.shape, self._in_offset,
                                    self.stops[-1]))
        # Check to see if we can process the next chunk and dump outputs
        while self._idx < len(self.starts) and \
                self._in_offset >= self.stops[self._idx]:
            start, stop = self.starts[self._idx], self.stops[self._idx]
            this_len = stop - start
            this_window = self._window.copy()
            if self._idx == len(self.starts) - 1:
                this_window = np.pad(
                    self._window, (0, this_len - len(this_window)), 'constant')
                for offset in range(self._step, len(this_window), self._step):
                    n_use = len(this_window) - offset
                    this_window[offset:] += self._window[:n_use]
            if self._idx == 0:
                for offset in range(self._n_samples - self._step, 0,
                                    -self._step):
                    this_window[:offset] += self._window[-offset:]
            logger.debug('    * Processing %d->%d' % (start, stop))
            this_proc = [in_[..., :this_len].copy()
                         for in_ in self._in_buffers]
            if not all(proc.shape[-1] == this_len == this_window.size
                       for proc in this_proc):
                raise RuntimeError('internal indexing error')
            outs = self._process(*this_proc, **kwargs)
            if self._out_buffers is None:
                max_len = np.max(self.stops - self.starts)
                self._out_buffers = [np.zeros(o.shape[:-1] + (max_len,),
                                              o.dtype) for o in outs]
            for oi, out in enumerate(outs):
                out *= this_window
                self._out_buffers[oi][..., :stop - start] += out
            self._idx += 1
            if self._idx < len(self.starts):
                next_start = self.starts[self._idx]
            else:
                next_start = self.stops[-1]
            delta = next_start - self.starts[self._idx - 1]
            for di in range(len(self._in_buffers)):
                self._in_buffers[di] = self._in_buffers[di][..., delta:]
            logger.debug('    - Shifting input/output buffers by %d samples'
                         % (delta,))
            self._store(*[o[..., :delta] for o in self._out_buffers])
            for ob in self._out_buffers:
                ob[..., :-delta] = ob[..., delta:]
                ob[..., -delta:] = 0.


def _check_cola(win, nperseg, step, window_name, tol=1e-10):
    """Check whether the Constant OverLap Add (COLA) constraint is met."""
    # adapted from SciPy
    binsums = np.sum([win[ii * step:(ii + 1) * step]
                      for ii in range(nperseg // step)], axis=0)
    if nperseg % step != 0:
        binsums[:nperseg % step] += win[-(nperseg % step):]
    const = np.median(binsums)
    deviation = np.max(np.abs(binsums - const))
    if deviation > tol:
        raise ValueError('segment length %d with step %d for %s window '
                         'type does not provide a constant output '
                         '(%g%% deviation)'
                         % (nperseg, step, window_name,
                            100 * deviation / const))
    return const


class _Storer(object):
    """Store data in chunks."""

    def __init__(self, *outs, picks=None):
        for oi, out in enumerate(outs):
            if not isinstance(out, np.ndarray) or out.ndim < 1:
                raise TypeError('outs[oi] must be >= 1D ndarray, got %s'
                                % (out,))
        self.outs = outs
        self.idx = 0
        self.picks = picks

    def __call__(self, *outs):
        if (len(outs) != len(self.outs) or
                not all(out.shape[-1] == outs[0].shape[-1] for out in outs)):
            raise ValueError('Bad outs')
        idx = (Ellipsis,)
        if self.picks is not None:
            idx += (self.picks,)
        stop = self.idx + outs[0].shape[-1]
        idx += (slice(self.idx, stop),)
        for o1, o2 in zip(self.outs, outs):
            o1[idx] = o2
        self.idx = stop


def _ensure_int(x, name='unknown', must_be='an int', *, extra=''):
    """Ensure a variable is an integer."""
    # This is preferred over numbers.Integral, see:
    # https://github.com/scipy/scipy/pull/7351#issuecomment-299713159
    extra = f' {extra}' if extra else extra
    try:
        # someone passing True/False is much more likely to be an error than
        # intentional usage
        if isinstance(x, bool):
            raise TypeError()
        x = int(operator.index(x))
    except TypeError:
        raise TypeError(f'{name} must be {must_be}{extra}, got {type(x)}')
    return x
