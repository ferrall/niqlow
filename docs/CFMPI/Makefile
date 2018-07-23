
PC := TRUE
ifeq ($(PC),TRUE)
	SED := "C:\Program Files (x86)\GNUWin32\bin\sed.exe"
else
	SED := sed
endif

tweak:
	$(SED) -i -f oxdoc.sed *.html
	erase sed*.
