CXX = g++
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
CXXFLAGS = -std=c++17 -Wall \
			-Iinclude \
			-I../../../UtilityCode/Include \
			$(ROOTCFLAGS)
TARGET = QuickDrawForest
SRC = 	$(TARGET).cpp \
		$(wildcard src/*.cpp) \
		$(wildcard ../../../UtilityCode/Source/*.cpp)
OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(ROOTLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f src/*.o $(TARGET)