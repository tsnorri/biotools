-include ../local.mk
include ../common.mk

OBJECTS	= cmdline.o skip_lines.o


all: skip_lines

clean:
	$(RM) $(OBJECTS) skip_lines

skip_lines: $(OBJECTS)
	$(CXX) $(COMMON_FLAGS) -o skip_lines $(OBJECTS)

skip_lines.cc: cmdline.c
