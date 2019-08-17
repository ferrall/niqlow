
PC := TRUE
ifeq ($(PC),TRUE)
	SED := "C:\Program Files (x86)\GNUWin32\bin\sed.exe"
	COPY := copy
    ERASE := erase
else
	SED := sed
	COPY := cp
    ERASE := rm
endif

tweak:
	$(SED) -i -f oxdoc.sed *.html
	$(SED) -i -f oxdoc1.sed DDP/*.html
	$(SED) -i -f oxdoc1.sed FiveO/*.html
	$(SED) -i -f oxdoc1.sed CFMPI/*.html
	$(SED) -i -f oxdoc1.sed About/*.html
	$(SED) -i -f oxdoc1.sed examples/*.html
	$(SED) -i -f oxdoc2.sed examples/DDP/*.html
	$(SED) -i -f oxdoc2.sed examples/FiveO/*.html
	$(SED) -i -f oxdoc2.sed examples/misc/*.html
	$(SED) -i -f oxdoc2.sed examples/replications/*.html
	$(ERASE) index.html
	$(COPY) niqlow.ox.html index.html
	$(ERASE) /S sed*.
