The HamiltonFastMarching can be used from Python in two ways:

1- using the library PythonHFM, compiled with pybind11 
2- using the executable FileHFM, which has no external dependencies

The two methods provide similar functionality. 
1 offers faster input/output, while 2 may be easier to compile.

When compiling with 1, you need to correctly set the path to the python executable.
Indeed, there are usually multiple python distributions on a system, and the PythonHFM library will not work if incorrectly linked. Type "which python" in a terminal to find about this path.