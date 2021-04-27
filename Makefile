#################################################################################
#		  NAME: MAKEFILE for nucleotide-search
# 		AUTHOR: David Rich
#################################################################################

# COMPILERS TOOLS:
SHELL 	:= /bin/sh
CC 		:= gcc
CXX 		:= g++
CP 		:= cp 

TARGET 	:= nucleotide-search-opt.out

SOURCE 	:= nucleotide-search.cpp  \
				Clock.hpp

CXX_FLAGs := -O3 -march=native -std=c++17

SOURCE_DIR := src
SOURCES   = $(SOURCE:%=$(SOURCE_DIR)/%)

$(TARGET): $(SOURCES)
	$(CXX) $(CXX_FLAGS) $(SOURCES) -o $(@) 
	
$(SOURCES): 
	$(CXX) $(@)