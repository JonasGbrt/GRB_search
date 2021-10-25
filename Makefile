CC := g++
CFLAGS += -g -Wall -O2 -I./include 
ROOT_FLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)
CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

.PHONY : all clean
all: GRB

GRB: GRB.o 
	$(CC) -o $@ $^ $(ROOT_LIBS) 


%.o: %.cpp
	$(CC) $(CFLAGS) $(ROOT_FLAGS) -c -o $@ $<

clean:
	@echo "Cleaning objects ..."
	@rm -f *.o GRB
	@rm -f *~

