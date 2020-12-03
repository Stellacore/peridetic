
# Generate API documentation with 'doxygen' tool
find_package(Doxygen)
if (DOXYGEN_FOUND)

	# Variables for expansion in Doxyfile template
	set(DoxyProjName "${PerideticProjName}")
	set(DoxyProjVersion "${PerideticProjVersion}")
	set( DoxyProjTagline
		"Easy Accurate Geodetic form/into Cartesian Coordinate Transformation"
		)

	# Default enable doxygen build option
	option(BUILD_DOC "Build API documentation" ON)

	# Edit template Doxyfile to insert project specific info
	set(DoxTemplate ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
	set(DoxOutFile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
	configure_file(${DoxTemplate} ${DoxOutFile} @ONLY)

	# command to invoke doxygen build
	add_custom_target(
		docs ALL
		COMMAND ${DOXYGEN_EXECUTABLE} ${DoxOutFile}
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating documentation with Doxygen utility"
		VERBATIM
		)

	install(
		DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
		DESTINATION
			${CMAKE_INSTALL_DOCDIR}/${local_LIBTARGET}
		)

else (DOXYGEN_FOUND)

	message("Doxygen tool required for building library API documentation")

endif (DOXYGEN_FOUND)

