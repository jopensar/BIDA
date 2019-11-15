CXX ?= g++
CFLAGS ?= -std=c++14 -Wall -O2 -march=native -I..
LDFLAGS ?=
SRCS := $(shell find . -name '*.cpp')
OBJS := $(SRCS:%.cpp=%.o)
DEPS := $(SRCS:%.cpp=%.d)

.PHONY: clean all FORCE

all: aps

aps: $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) -o aps $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -MMD -c $< -o $@

clean:
	rm -f aps $(OBJS) $(DEPS)

-include $(DEPS)
