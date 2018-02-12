#-------------------------------------------------------------------  
#           M I N G W 32-bit  T H R E E S   D L L   M A K E F I L E  
#-------------------------------------------------------------------  
  
all : compile
  
#------------------------------------------------  
#  GCC
#------------------------------------------------  
CC = gcc
INC = ../../include/C
CFLAGS = -Wall -O2 -m64 -I$(INC) -c 
LIBNAME = threes_64.dll
OXLIB = $(INC)/liboxwin.a
CURLIB = $(INC)/libcurldll.a
CPATH = .

%.o : $(CPATH)/%.c
	$(CC) $(CFLAGS) $<

OBJECTS = threes.o

compile : $(OBJECTS)
	$(CC) -L$(INC) -shared -o $(LIBNAME) $^ $(OXLIB) $(CURLIB) threes.def 
