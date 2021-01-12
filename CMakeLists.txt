cmake_minimum_required( VERSION 3.0 )

project( RhapsodyWTA LANGUAGES CXX)

# Excute conan outside or not.
option(CONAN_INSTALL_MANUALLY "Run Conan install command manually." OFF)
# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
endif()

if(NOT CONAN_INSTALL_MANUALLY)
    include(${CMAKE_SOURCE_DIR}/scripts/conan.cmake)
    conan_check(REQUIRED)
    conan_cmake_run(
        CONANFILE conanfile.py
        INSTALL_FOLDER ${CMAKE_BINARY_DIR}
        BUILD_TYPE ${CMAKE_BUILD_TYPE}
        BUILD missing
    )
endif()

set(CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR} ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory( src )