
PC := TRUE
ifeq ($(PC),TRUE)
	SED := "C:\Program Files (x86)\GNUWin32\bin\sed.exe"
	ERASE := erase
else
	SED := sed
endif

tweak:
	$(SED) -i -f oxdoc.sed *.html
	$(ERASE) sed*.
