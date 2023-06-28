import numpy as np
cimport numpy as np
from cython.view cimport array as cvarray
from typing import List

cdef class ArrayDict(dict, np.lib.mixins.NDArrayOperatorsMixin):
    cdef np.ndarray __array
    cdef tuple __all_keys

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for k, v in self.items():
            if isinstance(v, dict):
                self[k] = type(self)(**v)
            elif isinstance(v, np.ndarray):
                self[k] = type(self)(**{str(i): v[i] for i in range(len(v))})

    @property
    def array(self):
        if self.__array is None:
            self.__array = self.__array__()
        return self.__array

    @property
    def all_keys(self):
        if self.__all_keys is None:
            self.__all_keys = inner_all_keys(self)
        return self.__all_keys

    @property
    def shape(self):
        return self.array.shape

    def __array__(self):
        return inner_array(self)

    def __all_keys__(self):
        cdef list keys = []

        inner_all_keys(self, keys)

        return tuple(map(tuple, keys))

    def combine_dims(self, dims: tuple[int]):
        dims = sorted(dims)
        new_data = {}
        for key_tuple in self.all_keys:
            value = self
            for key in key_tuple:
                value = value[key]
            new_key_tuple = []
            new_key = ''
            for i, key in enumerate(key_tuple):
                if i in dims:
                    new_key += f'{key}-'
                else:
                    new_key_tuple.append(key)
            new_key_tuple.append(new_key[:-1])
            new_data[new_key_tuple[0]] = new_data.get(new_key_tuple[0], {})
            current_level = new_data[new_key_tuple[0]]
            for key in new_key_tuple[1:-2]:
                current_level[key] = current_level.get(key, {})
                current_level = current_level[key]
            current_level[new_key_tuple[-2]] = current_level.get(new_key_tuple[-2], {})
            current_level[new_key_tuple[-2]][new_key_tuple[-1]] = value
        return ArrayDict(**new_data)

cdef inner_array(data):
    if isinstance(data, dict):
        return concatenate_arrays([np.ndarray([inner_array(d)])
                               for d in data.values() if d is not False],
                              axis=0)
    else:
        return data

cdef inner_all_keys(data, list keys=None, int lvl=0):
    if keys is None:
        keys = []
    if isinstance(data, (int, float, str, bool)) or np.isscalar(data):
        return
    elif isinstance(data, dict):
        if len(keys) < lvl + 1:
            keys.append(list(data.keys()))
        else:
            keys[lvl] += [k for k in data.keys() if k not in keys[lvl]]
        for d in data.values():
            inner_all_keys(d, keys, lvl+1)
    elif isinstance(data, np.ndarray):
        rows = range(data.shape[0])
        if len(keys) < lvl+1:
            keys.append(list(rows))
        else:
            keys[lvl] += [k for k in rows if k not in keys[lvl]]
        if len(data.shape) > 1:
            inner_all_keys(data[0], keys, lvl+1)
    else:
        raise TypeError(f"Unexpected data type: {type(data)}")
    return tuple(map(tuple, keys))


def concatenate_arrays(arrays: list[np.ndarray], axis: int = None) -> np.ndarray:
    if axis is None:
        axis = 0
        arrays = [np.expand_dims(ar, axis) for ar in arrays]

    while axis < 0:
        axis += max(ar.ndim for ar in arrays)

    max_shape = np.max(get_homogeneous_shapes(arrays), axis=0)

    modified_arrays = []
    for arr in arrays:
        if len(arr) == 0:
            continue
        arr_shape = list(max_shape)
        arr_shape[axis] = arr.shape[axis]

        nan_array = np.full(arr_shape, np.nan)

        indexing = [slice(None)] * arr.ndim
        for ax in range(arr.ndim):
            if ax == axis:
                continue
            indexing[ax] = slice(0, arr.shape[ax])
        nan_array[tuple(indexing)] = arr

        modified_arrays.append(nan_array)

    result = np.concatenate(modified_arrays, axis=axis)

    return result


def get_homogeneous_shapes(arrays):
    # Determine the maximum number of dimensions among the input arrays
    max_dims = max([arr.ndim for arr in arrays])

    # Create a list to store the shapes with a homogeneous number of dimensions
    homogeneous_shapes = []

    # Iterate over the arrays
    for arr in arrays:
        # Get the shape of the array
        # Handle the case of an empty array
        if len(arr) == 0:
            shape = (0,)
            dims = 1
        else:
            shape = arr.shape
            dims = arr.ndim

        # Pad the shape tuple with additional dimensions if necessary
        num_dims_to_pad = max_dims - dims
        shape += (1,) * num_dims_to_pad

        # Add the shape to the list
        homogeneous_shapes.append(shape)

    return homogeneous_shapes