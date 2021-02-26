#! /bin/bash
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


echo "\nRepo/Build Info:"

#
# **CHANGE** these to reflect local system
#
PerideticInstallDir=/tmpLocal/  # Here Peridetic installed in custom location
	# If Peridetic installed in standard system locations
	# then can omit this and associated ...PREFIX_PATH below.
PerideticRepoDir=/repos/periUsage/  # Path to consuming project (this) repo

echo "\nBuild:"

# Assuming Peridetic is installed on system (Here in /tmpLocal)
# **NOTE**
#     : Note use of an optimized built option to ensure Peridetic performance
cmake \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_PREFIX_PATH=${PerideticInstallDir=} \
	${PerideticRepoDir} \
	#

# build example consuming project
make -j `nproc`

echo "\nTarget Libraries:"

# show library dependencies for the example executable programs
ldd hello*

echo "\nTarget Sizes:"

# show program size
uname -s -r -m -o # environment info
stat --format="FileSize: %n %s [bytes]" hello*

echo "\nTarget Exec:"

# execute (e.g. confirm built and that there's no optimization/removal)
./helloSansPeri
./helloWithPeri

