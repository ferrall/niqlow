Version := OxMetrics7
PC := TRUE
ifeq ($(PC),TRUE)
	OX := "C:\Program Files (x86)\$(Version)\Ox\bin\oxl.exe"
	OXFLAGS := -DOX7
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

vpath %.ox .
vpath %.h $(INC)
vpath %.oxo $(INC)
vpath %.ox.html ../../doc/Shared

Shared.oxo  :

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	$(COPY) $@ $(INC)
	$(ERASE) $@

.PHONY : clean
clean:
    $(info $(CURDIR))
	$(ERASE) *.oxo
	$(ERASE) $(INC)\Shared.oxo

.PHONY : document
document:
	$(ERASE) ..\..\doc\Shared\*.ox.html
	$(OXDOC) -include $(INC) *.ox

.PHONY : install
install:
	$(COPY) *.oxo $(INC)
