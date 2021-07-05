
add_subdirectory(yocto)
add_subdirectory(src)
include_directories(nanosvg/src)
add_subdirectory(realtime)
#add_subdirectory(eigen)

if(YOCTO_OPENGL)
add_subdirectory(yocto_gui)
endif(YOCTO_OPENGL)
