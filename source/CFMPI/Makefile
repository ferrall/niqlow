
CC = mpicc
INC := ../../include
DOC := ../../doc/CFMPI
OXDEV := /home/ferrallc/bin/OxMetrics7/ox/include
C64PATH=-L$(MPI)/lib -I$(MPI)/include -I$(INC) -I$(OXDEV)
CFLAGS = -fPIC -Wall -O2  -m64 -D__cdecl= $(C64PATH) -lmpi -c
SHFLAGS = -m64 -shared

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
	copy $(INC)/CFMPI.ox .
	copy $(INC)/MPIinterface.ox .
	copy $(INC)/useMPI.ox .
	ifeq ($(OS), Windows)
		"C:\Program Files (x86)\oxdoc\bin\oxdoc.bat" -include $(INC) -uplevel CFMPI.ox MPIinterface.ox InstallAndUse MPI_FAQ
		erase CFMPI.ox MPIinterface.ox useMPI.ox
	endif
	${MAKE} -C $(DOC) tweak
