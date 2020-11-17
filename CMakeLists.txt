# The name of our project
project(smash-analysis NONE C)

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
   message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# Minimum cmake version this is tested on
cmake_minimum_required(VERSION 3.5.1)

## Tell cmake where to find our modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(FindPythonModules)

find_package(PythonInterp 2.7 REQUIRED)

find_package(PythonLibs 2.7)
find_package(Cython)
if (PYTHONLIBS_FOUND AND CYTHON_FOUND)
    message(STATUS "Found PythonLibs and Cython, therefore can boost some code with cython.")
endif()

find_python_module(numpy VERSION 1.14.0 REQUIRED)
find_python_module(scipy VERSION 0.10.1 REQUIRED)
find_python_module(matplotlib VERSION 1.3.1 REQUIRED)
find_python_module(argparse VERSION 1.1 REQUIRED)
#find_python_module(yaml VERSION 3.11 REQUIRED)
#find_python_module(pandas VERSION 0.14.1 REQUIRED)

find_program(SMASH smash PATHS ${SMASH_PATH})
if(NOT SMASH)
    message(FATAL_ERROR
        "SMASH not found. Please specify a path to the SMASH build directory "
        "by passing '-DSMASH_PATH=...' to cmake.")
else()
    message(STATUS "Found SMASH: ${SMASH}")
    # copy executable
    file(COPY ${SMASH} DESTINATION ${CMAKE_BINARY_DIR})
endif()

# generate git statistics, requires gitstats in the path
find_program(GITSTATS gitstats)
if(GITSTATS)
   message(STATUS "gitstats found at ${GITSTATS}")
   add_custom_target(stats
                     COMMAND ${GITSTATS} ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/gitstats
                     COMMENT "Generating git statistics")
endif()


# Decide whether to plot experimental data as well. Note that in this case,
# the path to the experimental_data directory needs to be passed
if(EXP_DATA)
   option(WITH_EXP_DATA "Plot experimental data points as well." ON)
endif(EXP_DATA)

if(WITH_EXP_DATA)
  if(EXP_DATA)
    message(STATUS "Found experimental data: ${EXP_DATA}")
    message(STATUS "Comparison to experimental measurements will be plotted.")
  else()
    message(WARNING "Experimental data not found and will not be plotted. Provide path to experimental_data directory via
                     cmake .. -DEXP_DATA=path/to/exp/data")
  endif(EXP_DATA)
endif(WITH_EXP_DATA)

# 1) elementary tests
add_subdirectory(test/cross_sections)
add_subdirectory(test/angular_distributions)

# 2) box tests
add_subdirectory(test/detailed_balance)
add_subdirectory(test/elastic_box)
add_subdirectory(test/densities)

# 3) collider tests
add_subdirectory(test/FOPI_pions)
add_subdirectory(test/dileptons)
add_subdirectory(test/energy_scan)
