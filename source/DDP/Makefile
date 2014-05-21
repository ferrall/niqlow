Version := OxMetrics7
PC := TRUE
ifeq ($(PC),TRUE)
	OX := "C:\Program Files (x86)\$(Version)\Ox\bin\oxl.exe"
	OXFLAGS := -DOX7
	COPY := copy
	OXDOC := %OxDocRun%
	OXDOC := C:\Program Files (x86)\oxdoc\bin\oxdoc.bat
	SED := "C:\Program Files (x86)\GnuWin32\bin\sed.exe"
	ERASE := erase
	INC := ..\..\include
	DOC := ../../doc/DDP
else
	OX := /home/ferrallc/bin/$(Version)/ox/bin64/oxl
	OXFLAGS := -DMPI  -DOX7
	COPY := cp
	ERASE := rm
	INC := ../../include
	DOC := ../../doc/DDP
endif

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
	$(ERASE) $(INC)\DDP.oxo
	$(ERASE) $(INC)\DP.oxo
	$(ERASE) $(INC)\DPAuxiliary.oxo
	$(ERASE) $(INC)\StateVariable.oxo
	$(ERASE) $(INC)\Variables.oxo
	$(ERASE) $(INC)\Clock.oxo
	$(ERASE) $(INC)\Actions.oxo
	$(ERASE) $(INC)\TimeInvariant.oxo
	$(ERASE) $(INC)\Bellman.oxo
	$(ERASE) $(INC)\Methods.oxo
	$(ERASE) $(INC)\Outcomes.oxo
	$(ERASE) $(INC)\DPDebug.oxo
		
.PHONY : document
document:
	$(OXDOC) -include $(INC) -uplevel $(xtradoc) ..\Shared\Shared.ox
	${MAKE} -C $(DOC) tweak	

.PHONY : install
install:
	$(COPY) *.oxo $(INC)
