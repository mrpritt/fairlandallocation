override CFLAGS += -Wall -Wextra -pedantic -std=c++11 -MD -MP -O2
CPP = g++ 
LFLAGS =
INCLUDES =
LIBS = -static -lpthread -lstdc++ -lboost_program_options
SRC = $(wildcard *.cpp fmt/*.cpp lodepng/*.cpp)
HEADERS = $(wildcard *.h *.inl)
OBJS = $(SRC:.cpp=.o)
TARGET = proterra

all: $(TARGET)

debug: CFLAGS += -DHARD_DEBUG
debug: $(TARGET) 

release: CFLAGS:=$(filter-out -O2,$(CFLAGS))
release: CFLAGS += -DNDEBUG -O3
release: $(TARGET)
	
$(TARGET): $(OBJS)
	$(CPP) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LFLAGS) $(LIBS)

-include $(SRC:.cpp=.d)

%.o: %.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f *.o *.d $(TARGET) fmt/*.o fmt/*.d lodepng/*.o lodepng/*.d *.stackdump