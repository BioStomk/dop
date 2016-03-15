CC = g++
CFLAGS = -g -Wall

PROG = comptool
SRCS = burrows-wheeler-transform.cpp induced-sorting.cpp chaining.cpp comptool.cpp 
OBJS = $(SRCS:%.cpp=%.o)
DEPS = $(SRCS:%.cpp=%.d)


all: $(PROG)

$(PROG): $(OBJS)
	$(CC) -g -o $@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) -c -MMD -MP $<

-include $(DEPS)

.PHONY: clean
clean:
	-rm -f $(PROG) $(OBJS) $(DEPS)
