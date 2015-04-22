include ../Makefile.inc

INC := ..\..\include
DOC := ..\..\doc\Hybrids

vpath %.ox .:../Shared
vpath %.h $(INC)
vpath %.oxo $(INC)
vpath %.ox.html ../../doc/Hybrids

Hybrids.oxo : ImaiJainChing.oxo DPSystems.oxo

oxobjects = DDP.oxo Shared.oxo FiveO.oxo
ImaiJainChing.oxo : $(oxobjects)
DPSystems.oxo : $(oxobjects)

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	$(COPY) $@ $(INC)
	$(ERASE) $@


.PHONY : document
document:
	$(OXDOC) -include $(INC) -uplevel *.ox ../Shared/Help 
	${MAKE} -C $(DOC) tweak

