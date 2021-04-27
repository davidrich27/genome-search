#################################################################################
#		  NAME: MAKEFILE for nucleotide-search
# 		AUTHOR: David Rich
#################################################################################

# COMPILERS TOOLS:
SHELL 	:= /bin/sh
CC 		:= gcc
CXX 		:= g++
CP 		:= cp 

TARGET 	:= nucleotide-search

SOURCES 	:= nucleotide-search.cpp  \
				Clock.cpp

SOURCE_DIR := src 
SOURCES    := $(SOURCES:%=$(SOURCE_DIR)/%)

$(TARGET): $(SOURCES)
	$(CXX) -O3 $(SOURCES) -o $(@) 
	
$(SOURCES): 
	$(CXX) $(@)