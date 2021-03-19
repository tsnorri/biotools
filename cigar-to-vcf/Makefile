-include ../local.mk
include ../common.mk

OBJECTS	= cmdline.o main.o


all: cigar_to_vcf

clean:
	$(RM) $(OBJECTS) cigar_to_vcf

cigar_to_vcf: $(OBJECTS)
	$(CXX) $(COMMON_FLAGS) -o cigar_to_vcf $(OBJECTS)

main.cc: cmdline.c
