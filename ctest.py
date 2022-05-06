import ctypes
import pathlib
import sys

if __name__ == "__main__":
    libname = pathlib.Path().absolute() / "libraries/cfunctions.so"
    c_lib = ctypes.CDLL(libname)
    print(c_lib.cwindows.restype)
    values =c_lib.cwindows(b"g", 3, 3)


    print(type(values))