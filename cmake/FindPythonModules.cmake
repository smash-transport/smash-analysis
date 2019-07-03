function(find_python_module module)
    string(TOUPPER ${module} module_upper)
    set(options REQUIRED)
    set(oneValueArgs VERSION)
    cmake_parse_arguments(FIND_PY_${module_upper} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    if(${FIND_PY_${module_upper}_REQUIRED})
        set(PY_${module}_REQUIRED TRUE CACHE STRING
            "Whether Python module ${module} is required")
    endif()
    if(NOT PY_${module_upper}_FOUND)
        # A module's location is usually a directory, but for binary modules
        # it's a .so file.
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
            "import re, ${module}; print re.compile('/__init__.py.*').sub('',${module}.__file__)"
            RESULT_VARIABLE status
            OUTPUT_VARIABLE location
            ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT status)
            set(PY_${module_upper}_FOUND ${location} CACHE STRING
                "Location of Python module ${module}")
            execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
                "import re, ${module}; print re.compile('/__init__.py.*').sub('',${module}.__version__)"
                OUTPUT_VARIABLE version
                ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
            set(PY_${module_upper}_VERSION ${version} CACHE STRING
                "Version of Python module ${module}")
            if(FIND_PY_${module_upper}_VERSION)
                if(${version} VERSION_LESS ${FIND_PY_${module_upper}_VERSION})
                    message(SEND_ERROR "Found Python module ${module} version ${version}, but at least version ${FIND_PY_${module_upper}_VERSION} required!")
                else()
                    if(NOT FIND_PY_${module_upper}_REQUIRED)
                        set(optionally "optionally ")
                    else()
                        set(optionally "")
                    endif()
                    message(STATUS "Found Python module ${module} version ${version} (version ${FIND_PY_${module_upper}_VERSION} ${optionally}required).")
                endif()
            else()
                if(FIND_PY_${module_upper}_REQUIRED)
                    set(necessity "required")
                else()
                    set(necessity "optional")
                endif()
                message(STATUS "Found ${necessity} Python module ${module}.")
            endif()
        else()
            if(FIND_PY_${module_upper}_REQUIRED)
                message(SEND_ERROR "Required python module ${module} NOT found.")
            else()
                message(STATUS "Optional python module ${module} NOT found.")
            endif()
        endif()
    endif()
endfunction()
