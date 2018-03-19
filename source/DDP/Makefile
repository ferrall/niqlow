include ../Makefile.inc

INC := ..\..\include
DOC := ../../doc/DDP

vpath %.ox .:../Shared
vpath %.h $(INC):$(MENU)\source
vpath %.oxo $(INC)
vpath %.ox.html ../../doc/DDP

varsrc :=  Shared. DDPShared. ActionVariable. StateVariable. Clock. TimeInvariant. AuxiliaryValues.
varobj = $(varsrc:.=.oxo)
dppsrc := DP. Bellman. Methods. Data.
dppobj = $(dppsrc:.=.oxo)
dpsrc := DDPShared. Variables.
dpobj = $(dpsrc:.=.oxo)

xtradoc = DynProg GetStarted DDP.ox Error_And_Warning_Messages ../Shared/Help $(dpsrc:.=.ox) $(dppsrc:.=.ox)
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
	$(ERASE) $(INC)/DDPShared.oxo
	$(ERASE) $(INC)/AuxiliaryValues.oxo
	$(ERASE) $(INC)/StateVariable.oxo
	$(ERASE) $(INC)/Variables.oxo
	$(ERASE) $(INC)/Clock.oxo
	$(ERASE) $(INC)/Actions.oxo
	$(ERASE) $(INC)/TimeInvariant.oxo
	$(ERASE) $(INC)/Bellman.oxo
	$(ERASE) $(INC)/Data.oxo
	$(ERASE) $(INC)/Prediction.oxo
		
.PHONY : document
document:
	$(OXDOC) -include $(INC) -uplevel $(xtradoc) ../Shared/Shared.ox
	${MAKE} -C $(DOC) tweak	

.PHONY : install
install:
	$(COPY) *.oxo $(INC)
