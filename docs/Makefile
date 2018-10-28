
PC := TRUE
ifeq ($(PC),TRUE)
	SED := "C:\Program Files (x86)\GNUWin32\bin\sed.exe"
else
	SED := sed
endif

#	$(SED) -i -f oxdoc1.sed Hybrids/*.html
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
	erase /S sed*.
