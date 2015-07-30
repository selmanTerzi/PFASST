set(examples_to_install ${examples_to_install})

message(STATUS "  parareal")
include_directories(
    ${3rdparty_INCLUDES}
    ${FFTW_INCLUDE_PATH}
    ${pfasst_INCLUDES}
)

set(parareal_examples
    parareal_classic_serial
)

set(parareal_mpi_examples)

if(${pfasst_WITH_MPI})
    set(parareal_mpi_examples
        parareal_classic
        parareal_hybrid_full
        parareal_hybrid_partial
    )
    set(all_parareal_examples ${parareal_examples} ${parareal_mpi_examples})
else()
    set(all_parareal_examples ${parareal_examples})
endif()

foreach(example ${all_parareal_examples})
    add_executable(${example} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.cpp)
    if(NOT FFTW_FOUND)
        add_dependencies(${example} fftw3)
    endif()
    if(${pfasst_NUM_DEPENDEND_TARGETS} GREATER 0)
        add_dependencies(${example} ${pfasst_DEPENDEND_TARGETS})
    endif()
    target_link_libraries(${example}
        ${3rdparty_DEPENDEND_LIBS}
        ${FFTW_LIBRARIES}
    )
    if(pfasst_INSTALL_EXAMPLES)
        install(TARGETS ${example} RUNTIME DESTINATION bin)
    endif()
endforeach(example)

set(examples_to_install ${examples_to_install} ${all_parareal_examples} PARENT_SCOPE)
