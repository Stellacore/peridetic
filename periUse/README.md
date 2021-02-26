
	MIT License
	Copyright (c) 2020 Stellacore Corporation.

# Mini-Project Example - Independent Project Using Peridetic

of Peridetic use in external CMake build
# Mini-Project example of Peridetic use in external CMake build

This directory contains a trivial, but perhaps useful, compilable example
of a stand-alone andd independent 'project' that uses Peridetic capabilities.

## Description

This periUse 'project' builds two "hello-world" style programs. It uses
CMake to compile and link assuming Peridetic is already installed.

One purpose of this mini-project is to assess how much code/data space
overhead is added to consuming code in utilizing Peridetic capabilities.

The two programs themselves do very little. They allocate and set
data elements (to meaningless values), transform the data, and print
results. One program does not use Peridetic in any way. The other utilizes
Peridetic transformation operations (in both directions).

Comparing the sizes of the two exectuable programs provides an indication
of how much extra executable code is added by including Peridetic. (ref.
output of commands invoked from the runExample.sh bash script).


## Content

This "periUse" example includes:

	* README.md -- this file

	* runExample.sh -- script that illustrates how to configure CMake
		and use it to automatically find Peridetic and build an
		external consuming application (i.e. this periUse example)

	* CMakeLists.txt -- Minimalist CMake configuration for incorporating
		external (already installed) Peridetic into a consuming app.

	* Source code -- (helloCommon.h, helloSansPeri.cpp, helloWithPeri.cpp)
		Two 'main' programs used to comparing executable program size
		differences when either using, or not using, Peridetic in a program.

A screen shot "code diff" between the two programs is captured in the
[screen shot image](./periCodeSizeMains.png).

The build both with/without Peridetic use (both transform directions) produces

For example, on one system, with gcc 10.3, the executable programs have
the following sizes (exclusive of standard system/C++ dynamic link libraries).

	Target Sizes:
	Linux 4.15.0-122-generic x86_64 GNU/Linux
	FileSize: helloSansPeri 17208 [bytes]
	FileSize: helloWithPeri 21856 [bytes]


