Version := OxMetrics7
PC := TRUE
ifeq ($(PC),TRUE)
	OX := "C:\Program Files (x86)\$(Version)\Ox\bin\oxl.exe"
	OXFLAGS := -DOX7
	COPY := copy
	OXDOC := "C:\Program Files (x86)\oxdoc\bin\oxdoc.bat"
	ERASE := erase
	INC := ..\include
	XMP := ..\examples
	DOC := ../doc
else
	OX := /home/ferrallc/bin/$(Version)/ox/bin64/oxl
	OXFLAGS := -DMPI  -DOX7
	COPY := cp
	ERASE := rm
	INC := ../include
endif


vpath %.ox .:Shared:DDP:FiveO:Hybrids:$(INC)
vpath %.h $(INC)
vpath %.oxo $(INC)
vpath %.ox.html ../doc

oxobjects = Shared.oxo DDP.oxo FiveO.oxo CFMPI.ox Hybrids.oxo

niqlow.oxo  : $(oxobjects)

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	$(COPY) $@ $(INC)
	$(ERASE) $@

.PHONY : clean
clean:
	$(info $(CURDIR))
	$(ERASE) *.oxo
	$(ERASE) $(INC)\*.oxo

.PHONY : document
document:
	$(ERASE) ..\doc\default.html
	$(OXDOC) -uplevel -include $(INC) "niqlow is" "About niqlow"
	${MAKE} -C $(DOC) tweak
	${MAKE} -C About document
	${MAKE} -C DDP document
	${MAKE} -C FiveO document
	${MAKE} -C Hybrids document
	${MAKE} -C CFMPI document
	${MAKE} -C $(XMP) document
	
.PHONY : install
install:
	$(ERASE) $(INC)\*.oxo
	${MAKE} -C DDP
	${MAKE} -C FiveO
	$(COPY) FiveO\ParallelObjective.ox $(INC)
	${MAKE} -C Hybrids
	${MAKE} -C .
	$(ERASE) $(XMP)\*.oxo
	${MAKE} -C $(XMP)
	$(ERASE) $(XMP)\PeerCommunicationTestA.oxo
	$(ERASE) $(XMP)\ClientServerTestA.oxo
