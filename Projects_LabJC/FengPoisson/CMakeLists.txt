file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(FengPoisson ${sources} ${headers})

add_executable(FengPoisson ${FengPoisson})

target_link_libraries(FengPoisson pz)

