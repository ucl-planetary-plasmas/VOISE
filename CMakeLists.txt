
cmake_minimum_required(VERSION 3.19)
cmake_policy(VERSION 3.5)

#  Version and release date

set(VOISE_MAJOR 1)
set(VOISE_MINOR 0)
set(VOISE_PATCH 0)
set(VOISE_VERSION ${VOISE_MAJOR}.${VOISE_MINOR}.${VOISE_PATCH})


#   Set the project

project (VOISE VERSION ${VOISE_VERSION} LANGUAGES CXX)

#   Packaging

set(CPACK_GENERATOR "TBZ2")
set(CPACK_PACKAGE_VERSION_MAJOR ${VOISE_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VOISE_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VOISE_PATCH})
set(CPACK_SOURCE_GENERATOR "TBZ2")
set(CPACK_SOURCE_PACKAGE_FILE_NAME matlab-${VOISE_VERSION})
set(CPACK_SOURCE_IGNORE_FILES ".*.swp$" ".*.swo$" "~$" ".bz2$" ".gz$")

# Configure generators for binary installers and source packages.
include(CPack)
#include(UnitTests)
include(CheckCSourceCompiles)
#include(CreatePkgConfig)
# Define GNU standard installation directories
include(GNUInstallDirs)
# Helpers functions for creating config files that can be included by
# other projects to find and use a package
include(CMakePackageConfigHelpers)

include_directories(BEFORE SYSTEM ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
#  set(CMAKE_BUILD_TYPE Debug)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-DMX_COMPAT_64 -DMATLAB_MEXCMD_RELEASE=R2018a -DUSE_MEX_CMD   -D_GNU_SOURCE -DMATLAB_MEX_FILE -Wall -fexceptions -fPIC -O3 -fno-omit-frame-pointer -pthread -fwrapv -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-DMX_COMPAT_64 -DMATLAB_MEXCMD_RELEASE=R2018a -DUSE_MEX_CMD   -D_GNU_SOURCE -DMATLAB_MEX_FILE -Wall -fexceptions -fPIC -O0 -g -fno-omit-frame-pointer -pthread -fwrapv -fopenmp")


find_package(Matlab)

matlab_add_mex(NAME addSeedToVDBatch SRC cpp/mex/addSeedToVDBatch.cpp 
                    cpp/aux-functions/proposition2.cpp
										cpp/mex/pushVD.cpp
                    cpp/mex/grabVD.cpp
										cpp/skizException.cpp
										cpp/addSeed.cpp
										cpp/pointInRegion.cpp
										cpp/NSStar.cpp
										cpp/vd.cpp
										cpp/vd.cpp
										)

matlab_add_mex(NAME sopToWop SRC cpp/mex/sopToWop.cpp
                    cpp/mex/grabVD.cpp
										cpp/vd.cpp
										cpp/getRegion.cpp
										)

matlab_add_mex(NAME getCentroidSeedBatch SRC cpp/mex/getCentroidSeedBatch.cpp
                    cpp/mex/grabVD.cpp
										cpp/vd.cpp
										cpp/getRegion.cpp
										cpp/mex/grabW.cpp
										cpp/getCentroid.cpp
										)

matlab_add_mex(NAME getVDOpBatch SRC cpp/mex/getVDOpBatch.cpp
                    cpp/mex/pushVD.cpp
                    cpp/mex/grabVD.cpp
										cpp/skizException.cpp
										cpp/vd.cpp
										cpp/getRegion.cpp
										cpp/mex/grabW.cpp
										cpp/aux-functions/metrics.cpp
										cpp/getOp.cpp
										)

matlab_add_mex(NAME removeSeedFromVDBatch SRC cpp/mex/removeSeedFromVDBatch.cpp
                    cpp/mex/pushVD.cpp
                    cpp/mex/grabVD.cpp
										cpp/skizException.cpp
										cpp/removeSeed.cpp
										cpp/pointInRegion.cpp
										cpp/vd.cpp
										cpp/getRegion.cpp
										)

matlab_add_mex(NAME addSeedToVD SRC cpp/mex/addSeedToVD.cpp
                    cpp/aux-functions/proposition2.cpp
                    cpp/mex/pushVD.cpp
                    cpp/mex/grabVD.cpp
										cpp/skizException.cpp
										cpp/addSeed.cpp
										cpp/pointInRegion.cpp
										cpp/NSStar.cpp
										cpp/vd.cpp
										cpp/getRegion.cpp
										)

matlab_add_mex(NAME removeSeedFromVD SRC cpp/mex/removeSeedFromVD.cpp
                    cpp/mex/pushVD.cpp
                    cpp/mex/grabVD.cpp
										cpp/skizException.cpp
										cpp/removeSeed.cpp
										cpp/pointInRegion.cpp
										cpp/vd.cpp
										cpp/getRegion.cpp
										)

#   Installation
file(GLOB SOURCES "*.cpp" "*.h*")
install(FILES ${SOURCES} DESTINATION .)
install(FILES CMakeLists.txt DESTINATION .)
install(FILES ChangeLog DESTINATION .)



