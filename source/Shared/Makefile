include ../Makefile.inc

INC := ../../include
DOC := ../../docs/Shared

vpath %.ox .
vpath %.h $(INC)
vpath %.oxo $(INC)
vpath %.ox.html $(DOC)

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
	$(OXDOC) -include $(INC)  -uplevel Help *.ox

