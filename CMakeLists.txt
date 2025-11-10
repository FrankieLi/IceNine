# Modern CMakeLists for IceNine
cmake_minimum_required(VERSION 3.15)
project(IceNine VERSION 1.0 LANGUAGES CXX)


#----------------------------
set( IceNine_FileList  "Src/main.cpp" )

set( IceNine_FileList ${IceNine_FileList}  "Src/BoundarySelectionStrategies.cpp" )
set( IceNine_FileList ${IceNine_FileList}  "Src/BreadthFirstReconstructor.cpp" )
set( IceNine_FileList ${IceNine_FileList}  "Src/ConfigFile.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/CostFunctions.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/Detector.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/DetectorFile.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/DiscreteSearch.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/Driver.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/ExperimentSetup.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/ForwardSimulation.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/ImageData.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/InitFilesIO.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/main.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/OptimizationConfig.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/OrientationSearch.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/Peak.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/ReconstructionSetup.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/ReconstructionStrategies.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/Reconstructor.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/Sample.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/SearchDetails.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/Simulation.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/Utilities.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/XDMClient.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/XDMParallel.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/XDMRaster.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/XDMServer.cpp" )
set( IceNine_FileList ${IceNine_FileList} "Src/PaintGrid.cpp" )


#----------------------------

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add compiler flags to suppress deprecated warnings from Boost
add_compile_options(-Wno-deprecated-declarations)
add_compile_definitions(BOOST_ALLOW_DEPRECATED_HEADERS)

# Set build type to Release if not specified (CRITICAL for performance!)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Find MPI
find_package(MPI REQUIRED)

# Find Boost
find_package(Boost 1.43 REQUIRED)

# Force CMake to use only Homebrew Eigen3 and ignore 3rdParty directory
list(APPEND CMAKE_IGNORE_PATH "/Users/sfli/Research/3rdParty")
set(CMAKE_PREFIX_PATH "/opt/homebrew" ${CMAKE_PREFIX_PATH})

# Find Eigen3 (version 5.0 is backward compatible with 3.0)
find_package(Eigen3 REQUIRED NO_MODULE PATHS /opt/homebrew/share/eigen3/cmake NO_DEFAULT_PATH)

# Add XDM++ subdirectory
add_subdirectory(XDM++/libXDM)

# Create IceNine executable
add_executable(IceNine ${IceNine_FileList})

# Set include directories
target_include_directories(IceNine PRIVATE
    ${PROJECT_SOURCE_DIR}/XDM++/libXDM
    ${Boost_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_DIRS}
)

# Link Eigen3 (header-only library)
target_link_libraries(IceNine PRIVATE Eigen3::Eigen)

# Link libraries
target_link_libraries(IceNine PRIVATE
    XDM++
    ${MPI_CXX_LIBRARIES}
    ${Boost_LIBRARIES}
)

# Set MPI compile flags
target_compile_options(IceNine PRIVATE ${MPI_CXX_COMPILE_FLAGS}) 

