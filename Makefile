CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2
TARGET = ga_example
SOURCES = main2.cc problems.cc

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES)

clean:
	rm -f $(TARGET)

.PHONY:	clean
