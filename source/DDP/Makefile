include ../Makefile.inc

INC := ..\..\include
DOC := ../../doc/DDP


vpath %.ox .:../Shared
vpath %.h $(INC)
vpath %.oxo $(INC)
vpath %.ox.html ../../doc/DDP

varsrc := Shared. ActionVariable. StateVariable. Clock. TimeInvariant.
varobj = $(varsrc:.=.oxo)
dppsrc := DP. Bellman. Outcomes. Methods. DPDebug.
dppobj = $(dppsrc:.=.oxo)
dpsrc := DPAuxiliary. Variables.
dpobj = $(dpsrc:.=.oxo)

xtradoc = DDP.ox GetStarted GetStartedData $(varsrc:.=.ox) $(dppsrc:.=.ox) $(dpsrc:.=.ox)
xtrahtml = $(xtradoc:.ox=.ox.html)

DDP.oxo : $(dppobj)
DP.oxo : $(dpobj)
Variables.oxo : $(varobj)

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	$(COPY) $@ $(INC)
	$(ERASE) $@

.PHONY : clean
clean:
	$(info $(CURDIR))
	$(ERASE) $(INC)/DDP.oxo
	$(ERASE) $(INC)/DP.oxo
	$(ERASE) $(INC)/DPAuxiliary.oxo
	$(ERASE) $(INC)/StateVariable.oxo
	$(ERASE) $(INC)/Variables.oxo
	$(ERASE) $(INC)/Clock.oxo
	$(ERASE) $(INC)/Actions.oxo
	$(ERASE) $(INC)/TimeInvariant.oxo
	$(ERASE) $(INC)/Bellman.oxo
	$(ERASE) $(INC)/Methods.oxo
	$(ERASE) $(INC)/Outcomes.oxo
	$(ERASE) $(INC)/DPDebug.oxo
		
.PHONY : document
document:
	$(OXDOC) -include $(INC) -uplevel $(xtradoc) ../Shared/Shared.ox
	${MAKE} -C $(DOC) tweak	

.PHONY : install
install:
	$(COPY) *.oxo $(INC)
