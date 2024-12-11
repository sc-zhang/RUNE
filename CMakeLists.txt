cmake_minimum_required(VERSION 3.26)
project(RUNE)

set(CMAKE_CXX_STANDARD 11)

add_executable(rune src/main.cpp
        src/k_bin.cpp
        src/k_bin.h
        src/fasta_io.cpp
        src/fasta_io.h
        src/msg.h
        src/bin_io.cpp
        src/bin_io.h
        src/rune.h
        src/dumper.cpp
        src/dumper.h
        src/loader.cpp
        src/loader.h
        src/args.cpp
        src/args.h
)

find_package(ZLIB)
target_link_libraries(rune z)

install(TARGETS rune
        DESTINATION bin)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/scripts/rune_plot.py"
        DESTINATION scripts
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

