Version := OxMetrics7
OS := Windows

# list of subdirectories to make
DIRS := CFMPI DDP FiveO Hybrids Shared

ifeq ($(OS), Windows)
	OX := "C:\Program Files (x86)\$(Version)\Ox\bin\oxl.exe"
	OXFLAGS := -DOX7
	ECHO := echo
	CD := dir
	COPY := copy
	OXDOC := "C:\Program Files (x86)\oxdoc\bin\oxdoc.bat"
	ERASE := erase
	XMP := ..\examples
else ifeq ($(OS), Darwin)
	OX := oxl
	OXFLAGS := -DOX7
	ECHO := printf
	CD := cd
	COPY := cp
	ERASE := rm
else ifeq ($(OS), Linux)
	OX := /home/ferrallc/bin/$(Version)/ox/bin64/oxl
	OXFLAGS := -DMPI -DOX7
	ECHO := printf
	CD := cd
	COPY := cp
	ERASE := rm
endif
