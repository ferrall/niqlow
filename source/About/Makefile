
PC := TRUE
ifeq ($(PC),TRUE)
	OXDOC := "C:\Program Files (x86)\oxdoc\bin\oxdoc.bat"
	ERASE := erase
	INC := ..\..\include
	XMP := ..\examples
	DOC := ../../docs/About
else
	COPY := cp
	ERASE := rm
	INC := ../../include
endif
# ReleaseNotes
AboutDocs = Install License Citation Acknowledgements 

.PHONY : document
document:
	$(OXDOC) -uplevel -include $(INC) $(AboutDocs)
	
