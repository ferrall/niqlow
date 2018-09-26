include ../../source/Makefile.inc
DOC := ../../docs/examples/FiveO
INC := ../../include

vpath %.ox .
vpath %.h .:$(INC)
vpath %.oxo .:$(INC)
vpath %.ox.html $(DOC)

imports = AllTestOpt.
includes= GetStartedFiveO.ox 
exoxo = $(imports:.=.oxo)
exdoc =  $(includes) $(imports:.=.ox)

FiveOexamples.oxo : $(exoxo) $(includes)

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<	

.PHONY : document
document:
	$(OXDOC) -uplevel -include $(INC) $(exdoc)
