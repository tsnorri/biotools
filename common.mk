COMMON_FLAGS	?= -Wall -Werror -Wno-unused -O0 -g
GENGETOPT		?= gengetopt


%.o: %.cc
	$(CXX) -c $(COMMON_FLAGS) -std=c++17 -o $@ $<

%.o: %.c
	$(CC) -c $(COMMON_FLAGS) -std=c99 -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<"
