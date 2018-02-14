GENGETOPT	= gengetopt
CC			= clang-5.0
CXX			= clang++-5.0
CPPFLAGS	= -Wall -Werror -Wno-unused
CFLAGS		= -std=c99
CXXFLAGS	= -std=c++17
LDFLAGS		= -static-libstdc++ -static-libgcc

OBJECTS = cmdline.o process_imgt_gen_file.o

all: process_imgt_gen_file

clean:
	$(RM) $(OBJECTS)

process_imgt_gen_file: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS)

%.o : %.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

%.o : %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<"
