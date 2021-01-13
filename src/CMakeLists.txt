find_package(bamtools REQUIRED)
find_package(CLI11 REQUIRED)

set(HEADERS
    utils.h
)

set(SOURCES
    utils.cpp
    saturation.cpp
)

add_executable(saturation ${HEADERS} ${SOURCES})
target_include_directories(saturation PRIVATE ${bamtools_INCLUDES} ${CLI11_INCLUDES})
target_link_libraries(saturation ${bamtools_LIBRARIES} ${CLI11_LIBRARIES})