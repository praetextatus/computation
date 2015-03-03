CC = g++
SOURCES = $(wildcard *.cpp)
OBJS = $(SOURCES:.cpp=.o)
TARGET = tests
CFLAGS = -g -Wall -Wno-unknown-pragmas -std=c++0x
LFLAGS = -lboost_unit_test_framework -lboost_unit_test_framework-mt

$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(LFLAGS) -o $(TARGET) 


$(OBJS): %.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	touch $(OBJS)
	rm $(OBJS)
