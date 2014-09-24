file(GLOB SRCs_F src/*.f)

add_library(linpack ${SRCs_F})
