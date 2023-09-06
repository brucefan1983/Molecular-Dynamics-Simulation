CC = g++
CFLAGS = -std=c++14 -O3 -DDEBUG
INC = -I./

SOURCES = $(wildcard *.cpp)
ifdef OS # For Windows
OBJ = $(SOURCES:.cpp=.obj)
else
OBJ = $(SOURCES:.cpp=.o)
endif
HEADERS = $(wildcard *.h)

all: xiaomd
xiaomd: $(OBJ)
	$(CC) $(LDFLAGS) $^ -o $@ $(LIBS)
ifdef OS # for Windows
./%.obj: ./%.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@
else # for Linux
./%.o: ./%.cu $(HEADERS)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@
endif

clean:
ifdef OS
	del /s *.obj *.exp *.lib *.exe
else
	rm -f */*.o xiaomd
endif
