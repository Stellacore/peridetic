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

#
# === CMake info
#

cmake_minimum_required(VERSION 3.13)
set(PerideticProjName "peridetic")

set(PerideticVersion "0.1.0")

project(
	${PerideticProjName}
	VERSION ${PerideticVersion}
	LANGUAGES CXX
	)

# general configuration / dependencies
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# If any compile flags definitions - add here
add_definitions()


set(BUILD_FLAGS_FOR_CLANG
	""
	)
set(BUILD_FLAGS_FOR_GCC
	#
	-pthread
	-fomit-frame-pointer
	-fno-operator-names
	-fstrict-enums
	#
	-Wc++11-compat
	-Wall
	-Wextra
	-Wuninitialized
	-Winit-self
	-fno-nonansi-builtins
	-pedantic-errors
	)
set(BUILD_FLAGS_FOR_VISUAL
	""
	)


# target_link_libraries(... ${CMAKE_THREAD_LIBS_INIT})
#
# === Upstream CMake configurations
#


#
# === Configure project library
#

##@ set(peridetic_LIBVERSION ${PerideticProjName}-${PerideticVersion})

# === GNU standard paths
include(GNUInstallDirs)  # creates "CMAKE_INSTALL_<dir>"


# Project source code
add_subdirectory(include)  # public interface (and entire implementation)

# Specify documentation generation
add_subdirectory(doc)  # project documentation (via doxygen)

# Project test programs - via CTest
enable_testing()
add_subdirectory(eval)   # evaluation and experimentation programs
add_subdirectory(tests)  # quality and verification tests

# Example programs
add_subdirectory(examples)  # illustrative examples


#
# === Configure project unit tests
#

##@ add_subdirectory(test)


#
# === CMake configuration for downstream builds
#


#
# === Distribution packages
#

# Override default install location for testing this package
# In general default (no command here) will install to /usr/{lib,bin,...}
# - but here, set package to sandbox location for easier/safer testing
# set (DEB_PACKAGE_ROOT "/tmp/local/${PerideticProjName}")
# - default is "/usr" (e.g. for system packages)
# - specify explicitly as "/usr/local" as a developer local package
# - extract from command line
set (DEB_PACKAGE_ROOT "${CMAKE_INSTALL_PREFIX}")
set(CPACK_PACKAGING_INSTALL_PREFIX "${DEB_PACKAGE_ROOT}")

# list of packages to create (e.g. "DEB;RPM;...") if all tools installed
set(CPACK_GENERATOR "DEB")

#
# === Distribution package creation
#

	set(Peridetic_PACKAGE_VENDOR "https://github.com/Stellacore")
	set(Peridetic_PACKAGE_MAINTAINER "https://github.com/Stellacore/peridetic")

	# Specify generic cross-package distribuation properties
	# set(CPACK_PACKAGE_NAME ${PerideticProjName}) # same as default
	set(CPACK_PACKAGE_VENDOR "${Peridetic_PACKAGE_VENDOR}")
	set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
	set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
	set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
	set(CPACK_PACKAGE_VERSION_TWEAK ${PROJECT_VERSION_TWEAK})


	# == Debian package

		set(CPACK_PACKAGING_INSTALL_PREFIX "${DEB_PACKAGE_ROOT}")
		# Set package(Generator) specific characteristics
		set(CPACK_DEBIAN_PACKAGE_MAINTAINER "${Peridetic_PACKAGE_MAINTAINER}")
		# -- deps: none for peridetic
		#	set(APT_PACKAGE_DEPLIST
		#	"libeigen3-dev (>= ${EigenVersion}) \
		#	"
		#	)
		# set(CPACK_DEBIAN_PACKAGE_DEPENDS ${APT_PACKAGE_DEPLIST})


	# Must be included after the package variables are set
	include(CPack)

message("###")
message("### CMAKE_INSTALL_LIBDIR: " ${CMAKE_INSTALL_LIBDIR})
message("### CMAKE_INSTALL_BINDIR: " ${CMAKE_INSTALL_BINDIR})
message("### CMAKE_INSTALL_INCLUDEDIR: " ${CMAKE_INSTALL_INCLUDEDIR})
message("### CMAKE_INSTALL_DOCDIR: " ${CMAKE_INSTALL_DOCDIR})
message("###")
message("### CMAKE_INSTALL_FULL_LIBDIR: " ${CMAKE_INSTALL_FULL_LIBDIR})
message("### CMAKE_INSTALL_FULL_BINDIR: " ${CMAKE_INSTALL_FULL_BINDIR})
message("### CMAKE_INSTALL_FULL_INCLUDEDIR: " ${CMAKE_INSTALL_FULL_INCLUDEDIR})
message("### CMAKE_INSTALL_FULL_DOCDIR: " ${CMAKE_INSTALL_FULL_DOCDIR})
message("###")
message("### DEB_PACKAGE_ROOT: " ${DEB_PACKAGE_ROOT})
message("###")

