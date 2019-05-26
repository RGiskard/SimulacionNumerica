file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers} ${Commom_JC_GID} ${Commom_JC_Steklov} ${Adaptivity_JC_Refinement} ${Commom_JC_Head})

add_executable(TestAQP ${main_SRCS})

target_link_libraries(TestAQP pz)
