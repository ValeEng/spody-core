cmake_minimum_required(VERSION 3.10)
project(spody_core VERSION 1.0.0 LANGUAGES C)

# 1. Standard C (consigliato per precisione scientifica)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# 2. Elenco dei file sorgente (tutti in src/)
set(SOURCES
    src/spody_eclipse.c
    src/spody_ephemeris.c
    src/spody_harmonics.c
    src/spody_mapping.c
    src/spody_math.c
)

# 3. Creazione della libreria statica
# Produrrà libspody_core.a (Unix/macOS) o spody_core.lib (Windows)
add_library(spody_core STATIC ${SOURCES})

# 4. Configurazione dei percorsi degli header
# PUBLIC: la cartella include/ serve sia a chi compila la libreria
#         sia a chi la consuma (sia in-tree che dopo l'install)
target_include_directories(spody_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# 5. Link alla libreria matematica (essenziale per eclipse, harmonics, ecc.)
if(UNIX OR APPLE)
    target_link_libraries(spody_core PRIVATE m)
endif()

# 6. Regole di installazione
# Quando l'utente lancia `cmake --install . --prefix <path>`:
#   - la libreria finisce in <prefix>/lib
#   - tutti gli header in include/ finiscono in <prefix>/include
install(TARGETS spody_core
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

# Installa tutti gli header pubblici in un colpo solo
# (così se aggiungi nuovi header non serve toccare il CMakeLists)
install(DIRECTORY include/ DESTINATION include)
