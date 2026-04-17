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

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

set(SIM_SIPM_DIR "" CACHE PATH "Path to the SimSiPM directory")
if(NOT SIM_SIPM_DIR AND DEFINED ENV{SIM_SIPM_DIR})
  set(SIM_SIPM_DIR "$ENV{SIM_SIPM_DIR}" CACHE PATH "Path to the SimSiPM directory" FORCE)
endif()

find_path(PYTHIA8_INCLUDE_DIR NAMES Pythia8/Pythia.h)
find_library(PYTHIA8_LIBRARY NAMES pythia8)

if(NOT EXISTS "${SIM_SIPM_DIR}/include/SiPM.h")
  message(FATAL_ERROR "SiPM.h not found at ${SIM_SIPM_DIR}/include/SiPM.h")
endif()

if(NOT PYTHIA8_INCLUDE_DIR OR NOT PYTHIA8_LIBRARY)
  message(FATAL_ERROR "Pythia8 not found. Source the LCG view first.")
endif()

include(${Geant4_USE_FILE})

include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${ROOT_INCLUDE_DIRS}
  ${SIM_SIPM_DIR}/include
  ${PYTHIA8_INCLUDE_DIR}
)

link_directories(${SIM_SIPM_DIR}/build)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(HidraSim HidraSim.cc ${sources} ${headers})
target_link_libraries(HidraSim ${Geant4_LIBRARIES} sipm ${PYTHIA8_LIBRARY})

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

add_custom_target(G4HidraSim DEPENDS HidraSim)
install(TARGETS HidraSim DESTINATION bin)