
IF(BUILD_ASWIP)
    
    SET(ASWIP_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ASWIP}")

    SET(ASWIP_SOURCES src/sizes.F src/global.F src/global_3dvs.F src/boundaries.F 
                      src/mesh.F src/wind.F src/owiwind.F KDTREE2/kdtree2.F src/owi_ice.F 
                      wind/vortex.F wind/aswip_1.0.3.F )


    ADD_EXECUTABLE(aswip ${ASWIP_SOURCES})


    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        
        SET(ASWIP_COMPILER_FLAGS "${ASWIP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${ADDITIONAL_FLAGS_ASWIP}")
        
        SET_TARGET_PROPERTIES(aswip PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        
        TARGET_LINK_LIBRARIES(aswip netcdf netcdff )

    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        
        SET(ASWIP_COMPILER_FLAGS "${ASWIP_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG} ${ADDITIONAL_FLAGS_ASWIP}")
        SET(ASWIP_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(aswip PROPERTIES LINK_FLAGS ${ASWIP_LINKER_FLAGS} )
        
        TARGET_INCLUDE_DIRECTORIES(aswip PRIVATE ${CMAKE_SOURCE_DIR}/src)
        
        TARGET_LINK_LIBRARIES(aswip netcdf netcdff XdmfCore XdmfUtils Xdmf )

    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING) 

    TARGET_INCLUDE_DIRECTORIES(aswip PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_LINK_LIBRARIES(aswip version)

    SET_TARGET_PROPERTIES(aswip PROPERTIES COMPILE_FLAGS ${ASWIP_COMPILER_FLAGS})

    SET_TARGET_PROPERTIES(aswip PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/aswip_mod)

    INSTALL(TARGETS aswip RUNTIME DESTINATION bin)

ENDIF(BUILD_ASWIP)
