The HamiltonFastMarching can be used from Python in two ways:

1- using the library PythonHFM, compiled with boost-numpy (Experimental)
2- using the executable FileHFM, which has no external dependencies

The two methods provide similar functionality. 
1 offers faster input/output, while 2 may be easier to compile.

--- Compilation instructions for 1 (boost-numpy) ---
In CMake, one needs to select the link with the same Python library as Boost. 
See the discussion in http://stackoverflow.com/a/21052641

For me this is achieved by setting in CMake the following key:
PYTHON_LIBRARY = 
/usr/local/Cellar/python/2.7.13/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib

Important : Be careful to use the same version of python (2 or 3) for compiling the library and executing the scripts.
Failing to do results in the following error message (if compiling with python2 and executing scripts with python3)
ImportError: dynamic module does not define module export function (PyInit_PythonHFM_Isotropic2)

--------- Remarks on linking with python3 on Apple -------
By default, boost-python links with python2, and the python refers to python2, on this platform.
I have not tested python3 linking, but would be happy to now if anyone succeeds. 
Starting point:
brew uninstall --force boost-python
brew install boost-python --with-python3 --without-python