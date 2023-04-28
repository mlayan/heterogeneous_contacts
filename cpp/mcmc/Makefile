CC = g++
CFLAGS = -O2 -Wall -std=gnu++14
LDFLAGS =
INCLUDES =

SOURCES = main.cpp Household.cpp McmcObject.cpp utils.cpp
HEADERS = Household.h McmcObject.h utils.h
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = program.out

covid: $(EXECUTABLE) $(SOURCES) $(HEADERS)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(INCLUDES) $(OBJECTS) -o $@ -lstdc++fs

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean

clean:
	rm -rf *.o *.out
