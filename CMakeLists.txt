##**************************************************
## \file CMakeLists.txt
## \brief: CMakeLists of HidraSim project
## \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
## \start date: 7 July 2021
##**************************************************

# Setup the project
#
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(G4HidraSim)


# Find Geant4 package, activating all available UI and Vis drivers by default
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


# Define SIM_SIPM_DIR as a cache variable that can be configured
set(SIM_SIPM_DIR "" CACHE PATH "Path to the SimSiPM directory")

# Link SiPM sim
message(STATUS "SIM_SIPM_DIR from environment: $ENV{SIM_SIPM_DIR}")
include_directories(${SIM_SIPM_DIR}/include)
link_directories(${SIM_SIPM_DIR}/build)



# Setup Geant4 include directories and project include directories
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include
${ROOT_INCLUDE_DIRS})

# Locate project source and headers
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add the executable and link it to the Geant4 libraries
#
add_executable(HidraSim HidraSim.cc ${sources} ${headers})
#target_link_libraries(HidraSim ${Geant4_LIBRARIES})
#target_link_libraries(HidraSim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} sipm)
#target_link_libraries(HidraSim ${Geant4_LIBRARIES})
target_link_libraries(HidraSim ${Geant4_LIBRARIES} sipm)

# Copy all scripts to the build directory 
#
set(HidraSim_SCRIPTS
  HidraSim_gui.mac
  HidraSim_init_vis.mac
  HidraSim_run.mac
  )

foreach(_script ${HidraSim_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#For internal Geant4 use 
#
add_custom_target(G4HidraSim DEPENDS HidraSim)

# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS HidraSim DESTINATION bin)

##**************************************************
