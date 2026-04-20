cmake_minimum_required(VERSION 3.10)
project(spody_core VERSION 1.0.0 LANGUAGES C)

# 1. Standard C (opzionale, ma consigliato per precisione scientifica)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# 2. Elenco dei file sorgente (tutti in src/)
set(SOURCES
    src/eclipse.c
    src/ephemeris.c
    src/harmonics.c
    src/spody_mapping.c
)

# 3. Creazione della libreria
# Usiamo il nome 'spody_core' così produrrà libspody_core.a o spody_core.lib
add_library(spody_core STATIC ${SOURCES})

# 4. Configurazione dei percorsi degli Header
# PUBLIC: cartelle che servono anche a chi usa la libreria (include)
# PRIVATE: cartelle che servono solo a te per compilare i .c (src)
target_include_directories(spody_core PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
)

# 5. Link alla libreria matematica (essenziale per eclipse.c, harmonics.c, etc.)
if(UNIX OR APPLE)
    target_link_libraries(spody_core PRIVATE m)
endif()

# 6. REGOLE DI INSTALLAZIONE (Per creare il pacchetto da distribuire)
# Questo definisce dove vanno i file quando qualcuno lancia 'make install'
install(TARGETS spody_core DESTINATION lib)

# Installa l'header principale
install(FILES include/spody_core.h DESTINATION include)

# Installa gli header necessari che stanno in src
install(FILES 
    src/eclipse.h 
    src/ephemeris.h 
    src/harmonics.h 
    src/spody_mapping.h 
    DESTINATION include
)