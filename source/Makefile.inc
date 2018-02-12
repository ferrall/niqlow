Version := OxMetrics7
OS := Windows

# list of subdirectories to make
DIRS := CFMPI DDP FiveO Hybrids Shared C 

ifeq ($(OS), Windows)
#	OX := "C:\Program Files (x86)\$(Version)\Ox\bin\oxl.exe"
	OX := "C:\Program Files\$(Version)\Ox\bin64\oxl.exe"
	OXFLAGS := -DOX7 -v1 -d 
	ECHO := echo
	CD := dir
	COPY := copy
	OXDOC := "C:\Program Files (x86)\oxdoc\bin\oxdoc.bat"
	SED := "C:\Program Files (x86)\GnuWin32\bin\sed.exe"
	ERASE := erase
	XMP := ..\examples
else ifeq ($(OS), Darwin)
	OX := oxl
	OXFLAGS := -DOX7 -d 
	ECHO := printf
	CD := cd
	COPY := cp
	ERASE := rm
	SED := sed
else ifeq ($(OS), Linux)
	OX := /home/ferrallc/bin/$(Version)/ox/bin64/oxl
	OXFLAGS := -DMPI -DOX7 -d 
	ECHO := printf
	CD := cd
	COPY := cp
	ERASE := rm
	SED := sed
endif
