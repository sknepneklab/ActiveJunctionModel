# AJM (Active Junction Model)

This is the software package for simulating active junctions in the vertex model. The code accomanies
manuscript Generating active T1 transitions through mechanochemical feedback by Rastko Sknepnek, 
Ilyas Djafer-Cherif, Manli Chuai, Cornelis J. Weijer, and Silke Henkes submitted for consideration 
in eLife.

## Overview

The code is written in C++11. Analysis is does in Python 3 and MATLAB.

The main simulation code is written using modern C++ design. The user interact through a Python interface.
This is achieved by exposing key C++ to Python via the [pybind11](https://github.com/pybind/pybind11) library.

The software package produces output files in several formats (plain text, JSON, VTP, etc.). Notably, most of
the visualisation is done in [ParaView](https://www.paraview.org/), so having ParaView installed is
very helpful.

## Package Layout

src/         -- C++ source code
analysis/    -- Main analysis tools (Python)
extern/      -- External libraries 
make_conf/   -- Python scripts for generating various initial configuration
eLife/       -- Scripts and data used to produce figures in the manuscript

# Requirements

The code compiles on Linux, Mac OS X (tested up to OS X 10.15.7 Catalina), and in Windows Subsystem for Linux (WSL2) in Windows 11.

**Note**: In case there are issues with compiling the code on newer MacOS X systems, please use docker.

- Modern C++ compiler (tested with g++ and clang)
- CMake 3.x
- Boost (algorithm, lexical cast, property tree)
- VTK 7 or newer
- pybind11
- eigen 3 (optional)
- Python 3 (with NumPy, SciPy, shapely)


**Note**: One some Linux distributions and using VTK installed via anaconda, the code will not link do to a bug in the VTK package. 
In this case, it is recommended to compile and install VTK locally. In this case, one needs to set the CMake VTK_DIR variable to
/path/to/local/VTK/install/lib/cmake/vtk-9.1 (or other appropriate version).

# Compiling the code

Once all dependence are in place.

- mkdir build
- cd build
- ccmake ..
- make 

If the compilation is successful, a file called ajm.so (it may have different names on different systems) will be generated in the build directory.
This is the Python module. Please add the build directory to the PYTHONPATH environment variable or copy the .so file into the Python search path.

To test if everything works, in Python console write
```
import ajm
```

**Note**: In order for this work, the Python version has to be the same as the Python version against which the code was compiled.


