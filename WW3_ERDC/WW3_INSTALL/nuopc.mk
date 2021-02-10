#-----------------------------------------------
# NUOPC/ESMF self-describing build dependency
# makefile fragment for Wavewatch III
#-----------------------------------------------
# component module name
ESMF_DEP_FRONT := WMESMFMD
# component module path
ESMF_DEP_INCPATH := /home/xychen/ADC-WW3-Coupling-Dev/SrcCode/WW3/model/mod_MPI
# component module objects
ESMF_DEP_CMPL_OBJS := /home/xychen/ADC-WW3-Coupling-Dev/SrcCode/WW3/model/obj_MPI/libww3_multi_esmf.a(wmesmfmd.o)
# component object/archive list
ESMF_DEP_LINK_OBJS := /home/xychen/ADC-WW3-Coupling-Dev/SrcCode/WW3/model/obj_MPI/libww3_multi_esmf.a 
