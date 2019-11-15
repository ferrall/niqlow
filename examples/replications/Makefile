include ../../source/Makefile.inc
DOC := ../../docs/examples/replications
MDOC := ..\..\docs\examples\replications
INC := ../../include

vpath %.ox .
vpath %.h .:$(INC)
vpath %.oxo .:$(INC)
vpath %.ox.html $(DOC)

imports =RustEmet1987. RustEmet1987mle. WolpinJPE1984b. WolpinJPE1984. WolpinEmet1987. KennetJAE1994.  KeaneWolpinREStat1994. KeaneWolpinJPE1997. BerkovecSternEmet1991. GilleskieEmet1998. IJCEmet2009. AguirregabiriaMira2002.  Ahn1995. RosenzweigWolpinJPE1993. 
includes=RustEmet1987readdata.ox  
exoxo = $(imports:.=.oxo)
exdoc = $(examples:.=.ox)
exdoc =  $(includes) $(imports:.=.ox)

replications.oxo : $(exoxo) $(includes)

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<

.PHONY : document
document:
	$(OXDOC) -uplevel -include $(INC) $(exdoc)
	$(ERASE) /Q $(MDOC)\*.css
