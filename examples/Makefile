include ../source/Makefile.inc
DOC := ./
INC := ../include

vpath %.ox .
vpath %.h .:$(INC)
vpath %.oxo .:$(INC)
vpath %.ox.html .
#MortTest.   MVNormal. ReEmpExp.
#main.ox BerkovecSternEmet1991.ox KeaneWolpinREStat1994.ox TestSystemSolve.ox WolpinJPE1984b.ox WolpinJPE1984.ox WolpinEmet1987.ox RustEmet1987new.ox RustEmet1987readdata.ox RustEmet1987mle.ox KennetJAE1994.ox ReEmpExp.ox WstarTestb.ox GilleskieEmet1998.ox MNP.ox IJCEmet2009.ox AguirregabiriaMira2002.ox

examples =AllTest. AllTestOpt. RustEmet1987. RustEmet1987mle. WolpinJPE1984b. WolpinJPE1984. WolpinEmet1987. KennetJAE1994.  KeaneWolpinREStat1994. BerkovecSternEmet1991. WstarTestb. GilleskieEmet1998. MNP. IJCEmet2009. AguirregabiriaMira2002. TestGHK. StataMNP. ClientServerTestA. PeerCommunicationTestA.
exoxo = $(examples:.=.oxo)
exdoc = main.ox GetStarted.ox GetStartedData.ox RustEmet1987readdata.ox $(examples:.=.ox)

main.oxo : $(exoxo)

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	
.PHONY : clean
clean:
	$(info $(CURDIR))
	$(ERASE) /Q *.oxo

.PHONY : document
document:
	$(OXDOC) -uplevel -include $(INC) $(exdoc)
	$(SED) -i -f oxdoc.sed *.html						
	$(ERASE) sed*.

	
