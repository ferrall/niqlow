include ../Makefile.inc

INC := ..\..\include
DOC := ../../doc/FiveO

vpath %.ox .:../Shared
vpath %.h $(INC)
vpath %.oxo $(INC)
vpath %.ox.html ../../doc/FiveO

algsrc := Parameters. Algorithms. ParallelObjective.
algobj = $(algsrc:.=.oxo)
objsrc := Objective. BlackBox. Separable.
objobj = $(objsrc:.=.oxo)

FiveO.oxo : $(objobj)
Objective.oxo :  Shared.oxo $(algobj)
BlackBox.oxo : Objective.oxo
Separable.oxo : Objective.oxo

xtradoc = FiveO.ox GetStarted $(objsrc:.=.ox) $(algsrc:.=.ox) ../Shared/Shared.ox 

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	$(COPY) $@ $(INC)
	$(ERASE) $@
	
.PHONY : clean
clean:
	$(info $(CURDIR))
	$(ERASE) /Q *.oxo
	$(ERASE) /Q $(foreach f,$(objects) $(algobj),$(INC)\$f)

.PHONY : document
document:
	$(OXDOC) -include $(INC) -uplevel $(xtradoc)
	${MAKE} -C $(DOC) tweak

.PHONY : install
install:
	$(COPY) *.oxo $(INC)
