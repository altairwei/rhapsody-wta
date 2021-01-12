find_package(bamtools REQUIRED)

set(HEADERS
    utils.h
)

set(SOURCES
    utils.cpp
    seq_saturation.cpp
)

add_executable(seq-saturation ${HEADERS} ${SOURCES})
target_include_directories(seq-saturation PRIVATE ${bamtools_INCLUDES})
target_link_libraries(seq-saturation ${bamtools_LIBRARIES})