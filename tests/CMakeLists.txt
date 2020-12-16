#
#
# MIT License
#
# Copyright (c) 2020 Stellacore Corporation.
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


# Inventory of test programs
set(perideticTests
	
	testBuildHello # check build process
	testLocal # check development and testing tools
	testXforms # check fundamental transformation syntax/calling
	testCORS # check CORS data values parsing
	testAccuracy # check transformation external accuracy (vs CORS data)
	testMath # check various ellipsoid relationships

	)

# build exectutables
foreach (perideticTest ${perideticTests})

	# each test is independent executable program
	add_executable(${perideticTest} ${perideticTest}.cpp)

	# use may project build options
	target_compile_options(
		${perideticTest}
		PRIVATE
			$<$<CXX_COMPILER_ID:Clang>:${BUILD_FLAGS_FOR_CLANG}>
			$<$<CXX_COMPILER_ID:GNU>:${BUILD_FLAGS_FOR_GCC}>
			$<$<CXX_COMPILER_ID:MSVC>:${BUILD_FLAGS_FOR_VISUAL}>
		)

	# local and project include paths
	target_include_directories(
		${perideticTest}
		PRIVATE
			$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/
			$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../include/
        	$<INSTALL_INTERFACE:include/>
		)

	# dependency on project
	target_link_libraries(
		${perideticTest}
		PRIVATE
			peridetic::peridetic
		)

endforeach()

# add to "CTest" suite
foreach (perideticTest ${perideticTests})

	add_test("${perideticTest}"
		"${perideticTest}"
		)

endforeach ()


