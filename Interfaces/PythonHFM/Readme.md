# Calling the Hamiltonian Fast Marching library from Python

The HamiltonFastMarching code can be used from Python in two ways:

1. using the executable FileHFM, which has no external dependencies

2. using the library PythonHFM, compiled with pybind11

The two methods provide similar functionality.
Option 1. is easier to compile, whereas 2. offers faster input/output.

For option 1, change directory to ../FileHFM, and compile using CMake.

When compiling with 2, you need to correctly set the path to the python executable.
Indeed, there are usually multiple python distributions on a system, and the PythonHFM library will not work if incorrectly linked. Type "which python" in a terminal to find about this path.
