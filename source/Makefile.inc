#Version := OxMetrics7
#Version := OxMetrics8
#ExPath := Ox\bin64
Version := OxMetrics9
ExPath := Ox
OS := Windows

ifeq ($(OS), Windows)
	OX := "C:\Program Files\$(Version)\$(ExPath)\oxl.exe"
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
	XMP := ../examples
else ifeq ($(OS), Linux)
	OX := /home/ferrallc/bin/$(Version)/ox/$(ExPath)/oxl
	OXFLAGS := -DMPI -DOX7 -d
	ECHO := printf
	CD := cd
	COPY := cp
	ERASE := rm
	SED := sed
	XMP := ../examples
endif

# list of subdirectories to make Hybrids
DIRS := CFMPI Shared DDP FiveO C $(XMP)
