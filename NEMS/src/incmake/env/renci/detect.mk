########################################################################
#
# Main driver for supporting build environments in unknown clusters.
# Generally, this uses "uname" to detect the platform.
#
########################################################################

override uname_a=$(shell uname -a)

ifeq (,$(findstring croatan,$(uname_a)))
  NEMS_COMPILER?=intel
  $(call add_build_env, hatteras.$(NEMS_COMPILER),env/renci/hatteras.$(NEMS_COMPILER).mk)
endif
