set(examples_to_install ${examples_to_install})

message(STATUS "  parareal")
include_directories(
    ${3rdparty_INCLUDES}
    ${FFTW_INCLUDE_PATH}
    ${pfasst_INCLUDES}
)

set(parareal_examples
    parareal
)

set(parareal_mpi_examples)

if(${pfasst_WITH_MPI})
    set(parareal_mpi_examples
        parareal_mpi
    )
endif()

set(all_parareal_examples
    ${parareal_examples}
    ${parareal_mpi_examples}
)

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

if(${pfasst_WITH_MPI})
    include_directories(${MPI_CXX_INCLUDE_PATH})
    foreach(example ${parareal_mpi_examples})
        if(MPI_COMPILE_FLAGS)
            set_target_properties(${example} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
        endif()
        if(MPI_LINK_FLAGS)
            set_target_properties(${example} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
        endif()
        target_link_libraries(${example} ${MPI_CXX_LIBRARIES})
    endforeach(example)
endif()

set(examples_to_install ${examples_to_install} ${all_parareal_examples} PARENT_SCOPE)
