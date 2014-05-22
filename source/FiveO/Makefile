Version := OxMetrics7
PC := TRUE
ifeq ($(PC),TRUE)
	OX := "C:\Program Files (x86)\$(Version)\Ox\bin\oxl.exe"
	OXFLAGS := -DOX7
	COPY := copy
	OXDOC := C:\Program Files (x86)\oxdoc\bin\oxdoc.bat
	ERASE := erase
	INC := ..\..\include
	SED := "C:\Program Files (x86)\GNU\bin\sed.exe"
	DOC := ../../doc/FiveO
else
	OX := /home/ferrallc/bin/$(Version)/ox/bin64/oxl
	OXFLAGS := -DMPI  -DOX7
	COPY := cp
	ERASE := rm
	INC := ../../include
endif

vpath %.ox .:../Shared
vpath %.h $(INC)
vpath %.oxo $(INC)
vpath %.ox.html ../../doc/FiveO

partypes = Shared.oxo Parameters.oxo
algtypes = Shared.oxo Algorithms.oxo ParallelObjective.oxo
FiveO.oxo : Objective.oxo
Objective.oxo :  $(partypes) $(algtypes)

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	$(COPY) $@ $(INC)
	$(ERASE) $@
	
.PHONY : clean
clean:
	$(info $(CURDIR))
	$(ERASE) /Q *.oxo
	$(ERASE) /Q $(foreach f,$(objects) $(partypes) $(algtypes),$(INC)\$f)

.PHONY : document
document:
	$(OXDOC) -include $(INC) -uplevel GetStarted Parameters.ox Objective.ox ParallelObjective.ox Algorithms.ox ..\Shared\Shared.ox FiveO.ox
	${MAKE} -C $(DOC) tweak

.PHONY : install
install:
	$(COPY) *.oxo $(INC) 