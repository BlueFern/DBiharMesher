file(GLOB SRCs_F src/*.f)

add_library(blas ${SRCs_F})
