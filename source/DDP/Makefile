include ../Makefile.inc

INC := ..\..\include
DOC := ../../docs/DDP
MDOC := ..\..\docs\DDP


vpath %.ox .:../Shared
vpath %.h $(INC):$(MENU)\source
vpath %.oxo $(INC)
vpath %.ox.html $(DOC)

varsrc :=  Shared. DDPShared. ActionVariable. StateVariable. Clock. TimeInvariant. AuxiliaryValues.
varobj = $(varsrc:.=.oxo)
dppsrc := Methods. Bellman. Data.
dppobj = $(dppsrc:.=.oxo)
bellsrc := Variables. DDPShared. DP.
bellobj = $(bellsrc:.=.oxo)

xtradoc = GetStarted Guide DDP.ox DP.ox Variables.ox Bellman.ox TechAppendix Methods.ox Data.ox DDPShared.ox Error_And_Warning_Messages
xtrahtml = $(xtradoc:.ox=.ox.html)

DDP.oxo : $(dppobj)
Bellman.oxo : $(bellobj)
Variables.oxo : $(varobj)

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	$(COPY) $@ $(INC)
	$(ERASE) $@

.PHONY : document
document:
	$(OXDOC) -include $(INC) -uplevel $(xtradoc) ../Shared/Shared.ox
	$(COPY) $(MDOC)\default.html $(MDOC)\TOC.html
	$(ERASE) /Q $(MDOC)\*.css
