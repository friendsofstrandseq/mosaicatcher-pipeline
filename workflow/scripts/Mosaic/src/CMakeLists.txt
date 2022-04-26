cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project (Mosaicatcher CXX)


include(ExternalProject)


### Disable in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed. Please choose a different build path.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")



### Link Boost
find_package(Boost 1.50 COMPONENTS system date_time filesystem iostreams program_options REQUIRED)
include_directories (${Boost_INCLUDE_DIR})



### Link HTSlib using PackageConfig
find_package(PkgConfig REQUIRED) 

pkg_check_modules(HTSlib REQUIRED htslib) 
# set(HTSlib_LIBRARIES "${HTSlib_LDFLAGS}") 
# message(STATUS "   HTSlib include dirs: ${HTSlib_INCLUDE_DIRS}") 
# message(STATUS "   HTSlib libraries: ${HTSlib_LIBRARIES}")


# ### Set version

# # Get the current working branch
# execute_process(
#   COMMAND git rev-parse --abbrev-ref HEAD
#   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#   OUTPUT_VARIABLE GIT_BRANCH
#   OUTPUT_STRIP_TRAILING_WHITESPACE
# )

# # Get the latest abbreviated commit hash of the working branch
# execute_process(
#   COMMAND git log -1 --format=%h
#   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#   OUTPUT_VARIABLE GIT_COMMIT_HASH
#   OUTPUT_STRIP_TRAILING_WHITESPACE
# )


# # 1. Set Mosaicatcher version
# set (MOSAIC_VERSION_MAJOR 0)
# set (MOSAIC_VERSION_MINOR 3)

# # 2. Boost version from find_package (Boost_VERSION)

# # 3. Get HTSlib version from git
# execute_process(
#     COMMAND git describe --always
#     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/htslib
#     OUTPUT_VARIABLE HTSLIB_VERSION
#     OUTPUT_STRIP_TRAILING_WHITESPACE
# )
# set(HTSLIB_VERSION "1.3.1 ${HTSLIB_VERSION}")

# # 4. write version.hpp file
# configure_file (
#   "${PROJECT_SOURCE_DIR}/version.hpp.in"
#   "${PROJECT_BINARY_DIR}/version.hpp"
#   )

# # 5. Include build direcotry so version.hpp can be found
# include_directories("${PROJECT_BINARY_DIR}") # to find version.hpp



# ### Compiler flags
# add_definitions(-Wall -pedantic -Wno-unknown-pragmas)

 

# set(HEADER_FILES
#     count.hpp
#     counter.hpp
#     distribution.hpp
#     hmm.hpp
#     intervals.hpp
#     iocounts.hpp
#     program_options.hpp
#     segmentation.hpp
#     simulate.hpp
#     utils.hpp)

# ### Add executable
# add_executable (mosaic main.cpp ${HEADER_FILES})
# target_link_libraries (mosaic LINK_PUBLIC ${Boost_LIBRARIES})
# target_link_libraries (mosaic LINK_PUBLIC hts)

# ### Doxygen API documentation
# find_package(Doxygen)
# if(DOXYGEN_FOUND)
#     # write version into doxygen file
#     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
#     add_custom_target(doc ALL
#             COMMAND     ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
#             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
#             COMMENT "Generating API documentation with Doxygen" VERBATIM
#             )
# else()
#     message("Need Doxygen to produce API documentation. This is skipped.")
# endif(DOXYGEN_FOUND)
