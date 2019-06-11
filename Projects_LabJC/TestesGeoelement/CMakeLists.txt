file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})
add_executable(TestesGeoelement ${main_SRCS})

target_link_libraries(TestesGeoelement pz)

