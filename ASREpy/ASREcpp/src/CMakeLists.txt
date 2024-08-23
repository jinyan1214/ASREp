# Define minimum cmake version
cmake_minimum_required( VERSION 3.12 )

# Define project name
project(ASRElib VERSION 1.0)

# Declare an library
add_library(ASRElib SHARED ASRElib.cpp)

# Include the headers in inc
target_include_directories(ASRElib PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/inc")

# Create an INTERFACE target for the header-only library Eigen
add_library(Eigen INTERFACE)

# Specify the include directories for the header files of Eigen
target_include_directories(Eigen INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/inc/Eigen")

# Link Eigen library to ASRElib
target_link_libraries(ASRElib PRIVATE Eigen)


get_filename_component(parent_dir "${CMAKE_CURRENT_SOURCE_DIR}" DIRECTORY)

# If Windows add a MACRO _MSCV
IF(MSVC AND "${CMAKE_VS_PLATFORM_NAME}" MATCHES "(x64)")
    message("   ************************************")
    message("   Compiling on x64 Windows using MSVC")
    message("   ************************************")
    target_compile_definitions(ASRElib PRIVATE _MSCV)
    # Define the dir to install the compiled dll
    install( TARGETS ASRElib DESTINATION "${parent_dir}/bin/win32")
ELSEIF(MSVC)
    message("   ************************************")
    message("   Compiling on ${CMAKE_VS_PLATFORM_NAME} Windows using MSVC")
    message("   ************************************")
    target_compile_definitions(ASRElib PRIVATE _MSCV)
    install( TARGETS ASRElib DESTINATION "${parent_dir}/bin/win32")
ELSE(TRUE)
    message("   ************************************")
    message("   Compiling on unix-like")
    message("   ************************************")
ENDIF(MSVC AND "${CMAKE_VS_PLATFORM_NAME}" MATCHES "(x64)")



