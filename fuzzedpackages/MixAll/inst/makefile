#-----------------------------------------------------------------------
#     Copyright (C) 2004-2015  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
#
# Purpose:  Main makefile for the MixAll package
# Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
# Macros and definitions
RELATIVEMAKEVARS = makevars
include $(RELATIVEMAKEVARS)

#-----------------------------------------------------------------------
# main directories
MIXALL_PROJECTS_DIR  = projects/
MIXALL_INCLUDES_DIR  = include/
MIXALL_BIN_DIR = bin/
MIXALL_LIB_DIR = lib/
MIXALL_LIB = $(MIXALL_LIB_DIR)$(LIB_MIXALL_A)

# check if we are using symoblic paths
TOP_DIR      = ../../../
BIN_DIR      = $(TOP_DIR)$(MIXALL_BIN_DIR)
MAKEVARS     = $(TOP_DIR)$(RELATIVEMAKEVARS)
PROJECTS_DIR = $(TOP_DIR)$(MIXALL_PROJECTS_DIR)
INCLUDES_DIR = $(TOP_DIR)$(MIXALL_INCLUDES_DIR)

export RELATIVEMAKEVARS
export BIN_DIR
export MAKEVARS
export PROJECTS_DIR
export INCLUDES_DIR

#-----------------------------------------------------------------------
# Create the library. Build the objects of each projects and then add them
.PHONY: all lib obj
all: lib echo
lib: $(MIXALL_LIB)

# use to debug
echo:
	@echo $(BIN_DIR)
	@echo $(PROJECTS_DIR)
	@echo $(R_ARCH)
	@echo $(CXXFLAGS)
	@echo $(CXXPICFLAGS)
	@echo $(STK_CXXFLAGS)

#-----------------------------------------------------------------------
# projects paths
CLUSTERING_DIR= $(MIXALL_PROJECTS_DIR)Clustering/

SRC_CLUSTERING_DIR= $(CLUSTERING_DIR)src/

#-----------------------------------------------------------------------
# projects binaries
CLUSTERING_OBJS= $(patsubst $(SRC_CLUSTERING_DIR)%.cpp, $(MIXALL_BIN_DIR)%.o, $(wildcard $(SRC_CLUSTERING_DIR)*.cpp))

#-----------------------------------------------------------------------
# rules for compiling each projects
.PHONY: clustering

clustering:
	$(MAKE) -C $(SRC_CLUSTERING_DIR)

#-----------------------------------------------------------------------
# COMMENT or UNCOMMENT the projects you need for your package
# all binaries you need
obj: clustering

CPP_OBJS = \
  $(CLUSTERING_OBJS)

#-----------------------------------------------------------------------
# create the library adding the objects
$(MIXALL_LIB): obj
	$(AR) $(ARFLAGS) $@ $(CPP_OBJS)

#-----------------------------------------------------------------------
# rules for cleaning objects
.PHONY: clean clean_lib clean_obj clean_clustering

clean: clean_obj

clean_lib:
	$(RM) $(MIXALL_LIB)

clean_obj: \
  clean_clustering

clean_clustering:
	$(RM) $(CLUSTERING_OBJS)
