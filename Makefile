# define some Makefile variables for the compiler and compiler flags
# to use Makefile variables later in the Makefile: $()
#
#  -g          adds debugging information to the executable file
#  -Wall       turns on most, but not all, compiler warnings
#  -pedantic   turns on all warnings
#  -std=c++14  to use the latest C++ specifications

CC = g++
CFLAGS  = -g -Wall -pedantic -std=c++14
LIBS = -lboost_program_options -lm

# typing 'make' will invoke the first target entry in the file 
# (in this case the default target entry)
# you can name this target entry anything, but "default" or "all"
# are the most commonly used names by convention
#
default: motevo

# To create the executable file count we need the object files
# main.o parser.o and motevo.o
#
motevo:  parser.o motevo.o main.o
	$(CC) $(CFLAGS) $(LIBS) -o motevo.out parser.o motevo.o main.o

# To create the object file parser.o, we need the source
# files parser.cpp and parser.hpp
#
parser.o:  parser.cpp parser.hpp
	$(CC) $(CFLAGS) -c parser.cpp

# To create the object file motevo.o, we need the source files
# motevo.cpp, motevo.hpp and parser.hpp
#
motevo.o:  motevo.cpp motevo.hpp parser.hpp
	$(CC) $(CFLAGS) -c motevo.cpp

# To create the object file main.cpp we need the source files
# main.cpp motevo.hpp and parser.hpp
#
main.o:  main.cpp motevo.hpp parser.hpp
	$(CC) $(CFLAGS) -c main.cpp

# To start over from scratch, type 'make clean'.  This
# removes the executable file, as well as old .o object
# files and *~ backup files:
#
clean: 
	$(RM) count *.o *~
