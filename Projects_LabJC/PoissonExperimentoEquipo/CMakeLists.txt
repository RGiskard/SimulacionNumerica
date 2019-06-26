file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(PoissonEquipo ${sources} ${headers})

add_executable(PoissonEquipo ${PoissonEquipo})

target_link_libraries(PoissonEquipo pz)

