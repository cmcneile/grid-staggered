#
# Makefile to link in the grid library
# to some application code.
#
#  This makefile uses commands specific
#  to the gnu version of make.
#
#  This allows you to build a main program
#  against an installed grid library.
#
#  1) The location of the  program has to be set.
#  2) The main program you want to compile needs to be 
#     added to this makefile.
#
#

# The grid utility to extract the compile flags
gridinstall=/Users/cmcneile/code/grid/install
config=$(gridinstall)/bin/grid-config
#
#  flags from the grid utility program
#
CXXFLASGS= $(shell $(config) --cxxflags)   -I$(gridinstall)/include
LDFLAGS=$(shell $(config) --ldflags)
LIBS= $(shell $(config) --libs) -L$(gridinstall)/lib  -lGrid
COMPILER=$(shell $(config) --cxx)


cg.x: compute_staggered_meson.cc meson_routines.cc staggeredInclude.h
	$(COMPILER) -g -o $@ $(CXXFLASGS)  compute_staggered_meson.cc meson_routines.cc  $(LDFLAGS) $(LIBS) 

cg_block.x: compute_BLOCK_staggered_meson.cc meson_routines.cc staggeredInclude.h
	$(COMPILER) -g -o $@ $(CXXFLASGS)  compute_BLOCK_staggered_meson.cc meson_routines.cc  $(LDFLAGS) $(LIBS) 

io.x: dump_gauge_field.cc
	$(COMPILER) -g -o $@ $(CXXFLASGS) dump_gauge_field.cc   $(LDFLAGS) $(LIBS) 


mesonQED.x: compute_staggered_perturQED_meson.cc  meson_routines.cc staggeredInclude.h
	$(COMPILER) -g -o $@ $(CXXFLASGS) compute_staggered_perturQED_meson.cc  meson_routines.cc QED_meson_routines.cc  $(LDFLAGS) $(LIBS) 


