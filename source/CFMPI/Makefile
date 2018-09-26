include ../Makefile.inc

CC = mpicc
INC := ../../include
INCX := ..\..\include
DOC := ..\..\docs\CFMPI
DOCX := ../../docs/CFMPI
OXDEV := /home/ferrallc/bin/OxMetrics7/ox/include
C64PATH= -I$(INC) -I$(OXDEV) -I$(OXDEV)/dev $(LIB_DIRS)
CFLAGS = -fPIC -Wall -O2  -m64 $(C64PATH) -lmpi -c
SHFLAGS =  -shared

vpath %.ox $(INC)
vpath %.c .
vpath %.h $(INC)
vpath %.so $(INC)

CFMPI.so : CFMPI.o

%.o : %.c %.h
	$(CC)  $(CFLAGS) $<

%.so :
	$(CC) $(SHFLAGS) -o $@ $^
	mv $@ $(INC)

.PHONY : document
document:
	$(COPY)  $(INCX)\CFMPI.ox .
	$(COPY)  $(INCX)\MPIinterface.ox .
	$(COPY)  $(INCX)\useMPI.ox .
	$(OXDOC) -include $(INC) -uplevel CFMPI.ox MPIinterface.ox InstallAndUse MPI_FAQ GetStartedDIY
	$(ERASE) CFMPI.ox MPIinterface.ox useMPI.ox

