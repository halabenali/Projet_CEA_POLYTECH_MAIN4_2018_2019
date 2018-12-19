CXX = g++
FLAGS = -g -std=c++17
LIBFLAGS = -lblas -llapacke
SRC = main.cpp cholesky_facto.cpp file_API.cpp
OBJDIR = ./obj
OBJ = $(OBJDIR)/main.o $(OBJDIR)/cholesky_facto.o $(OBJDIR)/file_API.o
EXEC = time_exec_comp.exe

all: $(EXEC)

time_exec_comp.exe : $(OBJ)
	$(CXX) $(FLAGS) -o $@ $^ $(LIBFLAGS)

$(OBJDIR)/%.o: %.cpp
	$(CXX) $(FLAGS) -o $@ -c $<

.depend :
	g++ -MM $(SRC) > $@

-include .depend

clean:
	rm -f $(OBJ)
	rm -f $(EXEC)

.PHONY: clean all

