find_package(atlas REQUIRED)
find_package(eckit REQUIRED)

add_executable(atlasFVMNablaDriver atlasFVMNablaDriver.cpp)
target_link_libraries(atlasFVMNablaDriver atlas eckit DawnOptimizer DawnCompiler DawnUnittest GridTools::gridtools)
