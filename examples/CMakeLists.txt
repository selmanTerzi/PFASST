set(examples_to_install ${examples_to_install})

add_subdirectory(advection_diffusion)
add_subdirectory(scalar)
add_subdirectory(boris)
add_subdirectory(vanderpol)
add_subdirectory(parareal)

set(examples_to_install ${examples_to_install} PARENT_SCOPE)