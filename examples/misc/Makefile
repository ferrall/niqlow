include ../../source/Makefile.inc
DOC := ../../docs/examples/misc
MDOC := ..\..\docs\examples\misc
INC := ../../include

vpath %.ox .
vpath %.h .:$(INC)
vpath %.oxo .:$(INC)
vpath %.ox.html $(DOC)
#MortTest.   MVNormal. ReEmpExp. WstarTestb. DynamicWstar.

imports =MNP. 
includes=TestGHK.ox StataMNP.ox 
exoxo = $(examples:.=.oxo)
exdoc = $(includes) $(imports:.=.ox)

Miscexamples.oxo : $(exoxo) $(includes)

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	
.PHONY : document
document:
	$(OXDOC) -uplevel -include $(INC) $(exdoc)
	$(ERASE) /Q $(MDOC)\*.css
