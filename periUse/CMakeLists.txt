#
#
# MIT License
#
# Copyright (c) 2021 Stellacore Corporation.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
# KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
# AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
# IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
#

##
## Example for cmake build of an independent development effort that
## uses an installed version of Peridetic. E.g. assuming a Peridetic
## distro package has already been built and installed)
##



#
# === CMake info
#

# standard cmake items
cmake_minimum_required(VERSION 3.13)
set(PeriUseProjName "periUse") # arbitrary project name ...
set(PeriUseVersion "0.1.0")  # ... and it's arbitrary version ID

project(
	${PeriUseProjName}
	VERSION ${PeriUseVersion}
	LANGUAGES CXX
	)

# general configuration / dependencies
set(CMAKE_CXX_STANDARD 11)  # Peridetic header code is C++11 compliant
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# target_link_libraries
#
# === Upstream CMake configurations
#

# Perdidetic is used by (and is point of) example code
find_package(peridetic REQUIRED NO_MODULE)
message(Peridetic Found: ${peridetic_FOUND})
message(Peridetic Version: ${peridetic_VERSION})

set(targets

	helloWithPeri  # Includes peridetic.h and uses both transform directions
	helloSansPeri  # Independent of peridetic, but same code otherwise

	)

# Build example target programs
foreach(target ${targets})

	# standard cmake command for building a main program
	add_executable(${target} ${target}.cpp)

	# Here's how Peridetic is attached to development program
	# Incorporates the installed Peridetic header path into target build
	target_link_libraries(
		${target}
		PRIVATE
			peridetic::peridetic
		)

endforeach()

