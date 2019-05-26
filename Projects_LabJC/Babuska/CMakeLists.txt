file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Babuska ${sources} ${headers})

add_executable(Babuska ${Babuska})

target_link_libraries(Babuska pz)

