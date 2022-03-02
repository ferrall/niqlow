include ../Makefile.inc

INC := ..\..\include
DOC := ../../docs/FiveO
#:../Shared
vpath %.ox .
vpath %.h $(INC)
vpath %.oxo $(INC)
vpath %.ox.html $(DOC)

algsrc := Parameters. Algorithms. ParallelObjective.
algobj = $(algsrc:.=.oxo)
objsrc := Objective.
objobj = $(objsrc:.=.oxo)

FiveO.oxo : $(objobj)
#Shared.oxo 
Objective.oxo :  $(algobj)

xtradoc = FiveO.ox GetStarted $(objsrc:.=.ox) $(algsrc:.=.ox) Error_And_Warning_Messages ../Shared/Shared.ox

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
