-include ../local.mk
include ../common.mk

OBJECTS	= cmdline.o main.o


all: cigar_to_msa

clean:
	$(RM) $(OBJECTS) cigar_to_msa

cigar_to_msa: $(OBJECTS)
	$(CXX) $(COMMON_FLAGS) -o cigar_to_msa $(OBJECTS)

main.cc: cmdline.c
