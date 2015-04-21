CC = clang++
SOURCES = $(wildcard *.cpp)
OBJS = $(SOURCES:.cpp=.o)
TARGET = eigen
CFLAGS = -g -Wall -Wno-unknown-pragmas -std=c++0x

$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(LFLAGS) -o $(TARGET) 


$(OBJS): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	touch $(OBJS)
	rm $(OBJS)
