COMMON_FLAGS ?= -Wall -Werror -Wno-unused -O2 -g


%.o: %.cc
	$(CXX) -c $(COMMON_FLAGS) -std=c++14 -o $@ $<

%.o: %.c
	$(CC) -c $(COMMON_FLAGS) -std=c99 -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<"
