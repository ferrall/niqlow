include ../Makefile.inc

INC := ../../include

vpath %.ox .
vpath %.h $(INC)
vpath %.oxo $(INC)
vpath %.ox.html ../../doc/Shared

Shared.oxo  :

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	$(COPY) $@ $(INC)
	$(ERASE) $@

.PHONY : clean
clean:
	$(info $(CURDIR))
	$(ERASE) *.oxo
	$(ERASE) $(INC)\Shared.oxo

.PHONY : document
document:
	$(ERASE) ..\..\doc\Shared\*.ox.html
	$(OXDOC) -include $(INC) *.ox

.PHONY : install
install:
	$(COPY) *.oxo $(INC)
