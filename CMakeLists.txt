# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)


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

# This project is designed to be built outside the Insight source tree.
PROJECT(IceNine)

SET(CMAKE_CXX_COMPILER mpic++-mpich-gcc49)

include_directories("/Users/li31/LLNL/BES/Mesh_3/ExtLib/eigen")    # using eigen

include_directories("${PROJECT_SOURCE_DIR}/XDM++.strain/libXDM")    # define location of XDM++ (user download)
 add_subdirectory( XDM++.strain/libXDM )
SET(Boost_DIR /Users/Shared/Libraries/boost-trunk/boost)            # define location of boost 1.43 (user download)
include_directories(/opt/local/include/)


ADD_EXECUTABLE(IceNine ${IceNine_FileList} )
TARGET_LINK_LIBRARIES( IceNine XDM++  ) 

