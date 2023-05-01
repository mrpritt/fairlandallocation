override CFLAGS += -Wall -Wextra -std=c++14 -MD -MP -O3

LFLAGS = 
INCLUDES =
LIBS = -lpthread -lstdc++ -lboost_program_options
SRC = $(wildcard *.cpp ../../src/fmt/*.cpp)
HEADERS = $(wildcard *.h *.inl)
OBJS = $(SRC:.cpp=.o)
TARGET = generate

all: $(TARGET)

$(TARGET): $(OBJS)
	g++ $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LFLAGS) $(LIBS)

-include $(SRC:.cpp=.d)

%.o: %.cpp
	g++ $(CFLAGS) $(INCLUDES) -c $< -o $@

rebuild: clean all

clean:
	rm -f *.o *.d $(TARGET) ../../src/fmt/*.o ../../src/fmt/*.d 
