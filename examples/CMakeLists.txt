#
# Example build agains peridetic
#

#
# === Example program
#

set(perideticExamps

	lpaForXyz
	xyzForLpaRadians

	)

# build exectutables
foreach (perideticExamp ${perideticExamps})

	# each test is independent executable program
	add_executable(${perideticExamp} ${perideticExamp}.cpp)

	# use may project build options
	target_compile_options(
		${perideticExamp}
		PRIVATE
			$<$<CXX_COMPILER_ID:Clang>:${BUILD_FLAGS_FOR_CLANG}>
			$<$<CXX_COMPILER_ID:GNU>:${BUILD_FLAGS_FOR_GCC}>
			$<$<CXX_COMPILER_ID:MSVC>:${BUILD_FLAGS_FOR_VISUAL}>
		)

	# local and project include paths
	target_include_directories(
		${perideticExamp}
		PRIVATE
			$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/
			$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/../include/
        	$<INSTALL_INTERFACE:include/>
		)

	# dependency on project
	target_link_libraries(
		${perideticExamp}
		PRIVATE
			peridetic::peridetic
		)

endforeach()

