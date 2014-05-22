Version := OxMetrics7
PC := TRUE
ifeq ($(PC),TRUE)
	OX := "C:\Program Files (x86)\$(Version)\Ox\bin\oxl.exe"
	OXFLAGS :=  -DOX7
	COPY := copy
	OXDOC := C:\Program Files (x86)\oxdoc\bin\oxdoc.bat
	ERASE := erase
	INC := ..\..\include
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
	$(OXDOC) -include $(INC) -uplevel *.ox
