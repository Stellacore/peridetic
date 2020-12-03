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

# Create a header-only "library"
set(PerideticLibName "peridetic")

# Files comprising header installation
Set(PerideticHeaders

	peridetic.h   # public interface
	periDetail.h  # underlying implementation of peridetic.h

	)

# Define as a header-only interface
add_library(${PerideticLibName} INTERFACE)

# Alias for use in downstream CMakeLists
add_library(Peridetic::${PerideticProjName} ALIAS ${PerideticLibName})

# # add (header) sources
# target_sources(
# 	${PerideticLibName}
# 	INTERFACE
# 		${PerideticHeaders}
# 	)

# Export info
set(EXPORT_ALIAS_NAME ${PerideticProjName}::${PerideticLibName})
set(EXPORT_TARGETS_NAME ${PerideticProjName}Targets)
set(EXPORT_FILENAME ${EXPORT_TARGETS_NAME}.cmake)
add_library(
	${EXPORT_ALIAS_NAME}
	ALIAS
	${PerideticLibName}
	)

message("### CMAKE_MAJOR_VERSION: " ${CMAKE_MAJOR_VERSION})
message("### CMAKE_MINOR_VERSION: " ${CMAKE_MINOR_VERSION})
message("### CMAKE_PATCH_VERSION: " ${CMAKE_PATCH_VERSION})
message("### CMAKE_TWEAK_VERSION: " ${CMAKE_TWEAK_VERSION})
message("### PROJECT: " ${PROJECT})
message("### PROJECT_VERSION: " ${PROJECT_VERSION})
message("### PerideticLibName: " ${PerideticLibName})
message("### EXPORT_TARGETS_NAME: " ${EXPORT_TARGETS_NAME})
message("### EXPORT_FILENAME: " ${EXPORT_FILENAME})

			####

#[[
# install target (for use with downstream cmake find() command)
install(
	TARGETS ${PerideticLibName}
	EXPORT ${EXPORT_TARGETS_NAME}
	INCLUDES DESTINATION
		$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${PerideticProjName}>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	)
#]]

# install header files
install(
	FILES ${PerideticHeaders}
	DESTINATION
		${CMAKE_INSTALL_INCLUDEDIR}/${PerideticProjName}
	)

# CMAKE_INSTALL_INCLUDEDIR

			####

set(Peridetic_INCDIR include/${PerideticLibName})
install(
	TARGETS ${PerideticLibName}
	EXPORT ${EXPORT_TARGETS_NAME}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION
		$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${Peridetic_INCDIR}>
		$<INSTALL_INTERFACE:${Peridetic_INCDIR}>
	)

message("### Peridetic_INCDIR: " ${Peridetic_INCDIR})

install(
	EXPORT ${EXPORT_TARGETS_NAME}
	FILE ${EXPORT_FILENAME}
	NAMESPACE ${PerideticProjName}::
	DESTINATION
		${CMAKE_INSTALL_LIBDIR}/cmake/${PerideticLibName}
	)


##
## Creates Config.cmake and ConfigVersion.cmake files (that use install export)
##

include(CMakePackageConfigHelpers)

message("CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR})

configure_package_config_file(
	${CMAKE_SOURCE_DIR}/cmake/xConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/${PerideticProjName}Config.cmake
	INSTALL_DESTINATION
		${CMAKE_INSTALL_LIBDIR}/cmake/${PerideticLibName}
	)

write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/${PerideticProjName}ConfigVersion.cmake
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
)

install(
	FILES
		${CMAKE_CURRENT_BINARY_DIR}/${PerideticProjName}Config.cmake
		${CMAKE_CURRENT_BINARY_DIR}/${PerideticProjName}ConfigVersion.cmake
	DESTINATION
		${CMAKE_INSTALL_LIBDIR}/cmake/${PerideticLibName}
)

#[[
#]]

