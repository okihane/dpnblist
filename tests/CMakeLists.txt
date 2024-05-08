set(TEST_EXE dpnblist_test)

add_executable(${TEST_EXE}
    main.cpp
    test_box.cpp
    read_lmp.cpp
    test_cellSchAlgCPU.cpp
    test_cellSchAlgGPU.cpp
    test_octreeSchAlgCPU.cpp
    test_octreeSchAlgGPU.cpp
    #test_hashSchAlgCPU.cpp
    #test_hashSchAlgGPU.cpp
)

find_package(doctest QUIET)
if (NOT doctest_FOUND)
    set(DOCTEST_DIR ${CMAKE_SOURCE_DIR}/external/doctest-2.4.11)
    
endif()

set(HEADER_FILES
    src/nbList.h
    src/box.h
    src/nbSchAlg.h
    src/SchAlg/hashSchAlgGPU.cuh
    src/SchAlg/hashSchAlgCPU.h
    src/SchAlg/octreeSchAlgCPU.h
    src/SchAlg/octreeSchAlgGPU.cuh
    src/SchAlg/cellSchAlgCPU.h
    src/SchAlg/cellSchAlgGPU.cuh
)

target_include_directories(${TEST_EXE} PRIVATE
    ${DOCTEST_DIR}
    ${HEADER_FILES}
)

target_link_libraries(${TEST_EXE} PRIVATE nblist)
