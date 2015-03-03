CC = g++
SOURCES = $(wildcard *.cpp)
OBJS = $(SOURCES:.cpp=.o)
TARGET = solve2
CFLAGS = -g -Wall -Wno-unknown-pragmas -std=c++0x

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET) 

$(OBJS): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	touch $(OBJS)
	rm $(OBJS)
