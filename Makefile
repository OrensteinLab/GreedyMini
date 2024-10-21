CXX = g++
CXXFLAGS = -std=c++20 -O3 -march=native -I ~/boost_1_82_0
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
TARGET = my_program

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(OBJECTS) $(TARGET)
