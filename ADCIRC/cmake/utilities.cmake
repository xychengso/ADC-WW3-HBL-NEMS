
IF(BUILD_UTILITIES)
    SET(UTILITY_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${ADCIRC_OPTION_FLAGS}")
    
    ADD_EXECUTABLE(adccmp util/adccmp.F)
    ADD_EXECUTABLE(p15 wind/p15.F)
    ADD_EXECUTABLE(owi22 wind/owi22.F)
    ADD_EXECUTABLE(build13 util/build13.F)
    ADD_EXECUTABLE(buildstwave23 util/buildstwave23.F)
    ADD_EXECUTABLE(hot2asc util/hot2asc.F)
    ADD_EXECUTABLE(inflate util/inflate.F)
    ADD_EXECUTABLE(hstime util/hstime.F)
    
    SET_TARGET_PROPERTIES(adccmp PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(p15 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(owi22 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(build13 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(buildstwave23 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(hot2asc PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(inflate PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)
    SET_TARGET_PROPERTIES(hstime PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/util_mod)

    SET_TARGET_PROPERTIES(adccmp PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(p15 PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(owi22 PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(build13 PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(buildstwave23 PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(hot2asc PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(inflate PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})
    SET_TARGET_PROPERTIES(hstime PROPERTIES COMPILE_FLAGS ${UTILITY_COMPILER_FLAGS} ${ADDITIONAL_FLAGS_UTILITIES})

    INSTALL(TARGETS adccmp p15 owi22 build13 buildstwave23 hot2asc inflate hstime 
            RUNTIME DESTINATION bin)

ENDIF(BUILD_UTILITIES)
