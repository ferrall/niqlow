include ../source/Makefile.inc

vpath %.ox .
vpath %.oxh .
vpath %.oxo .

CFMPI.oxo : 

%.oxo : %.ox %.oxh 
	$(OX) $(OXFLAGS) -c $<
	$(ERASE) $@
