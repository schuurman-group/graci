"""
Module for facilitating passing of strings between
Python and compiled dlls
"""
import sys as sys
import ctypes as ctypes

def convert_ctypes(py_val, dtype=None):
    """convert a python array into a C-data type"""

    # note: the current approach is used based on:
    # https://bugs.python.org/issue27926
    # Namely, the default Python constructor is _very_ slow.
    # if that changes, then this function can change too...

    # there are fancier ways to do this by querying both the array and
    # machine environment, but this will do for now
    if dtype == 'int32':
        type_sym = 'i';i_size = 4; ctype_sym = ctypes.c_int32
    elif dtype == 'int64':
        type_sym = 'i';i_size = 8; ctype_sym = ctypes.c_int64
    elif dtype == 'double':
        type_sym = 'd';i_size = 8; ctype_sym = ctypes.c_double
    elif dtype == 'logical':
        type_sym = 'i';i_size = 4; ctype_sym = ctypes.c_bool
    elif dtype == 'string':
        ctype_sym = ctypes.c_char_p
    else:
        sys.exit('convert_ctypes does not recognize dtype='+str(dtype))

    if isinstance(py_val, str):
        return ctype_sym(py_val.encode('utf-8'))

    if isinstance(py_val, (float, int)):
        return ctype_sym(py_val)

    if py_val.size == 1:
        c_arr = (ctype_sym * py_val.size)(py_val)
    else:
        c_arr = (ctype_sym * py_val.size)(*py_val)

    return c_arr


