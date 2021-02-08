# ESMF self-describing build dependency makefile fragment

ESMF_DEP_FRONT     = adc_cap
ESMF_DEP_INCPATH   = /home/xychen/ADC-WW3-Coupling-Dev/SrcCode/ADCIRC/cpl/nuopc /home/xychen/ADC-WW3-Coupling-Dev/SrcCode/ADCIRC_INSTALL 
ESMF_DEP_CMPL_OBJS = 
ESMF_DEP_LINK_OBJS =  -L/home/xychen/ADC-WW3-Coupling-Dev/SrcCode/ADCIRC_INSTALL -ladc /home/xychen/ADC-WW3-Coupling-Dev/SrcCode/ADCIRC_INSTALL/libadc_cap.a  -L/home/xychen/ADC-WW3-Coupling-Dev/SrcCode/ADCIRC/work/  /home/xychen/ADC-WW3-Coupling-Dev/SrcCode/ADCIRC/work/libadc.a  
