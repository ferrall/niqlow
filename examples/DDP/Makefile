include ../../source/Makefile.inc
DOC := ../../docs/examples/DDP
INC := ../../include

vpath %.ox .
vpath %.h .:$(INC)
vpath %.oxo .:$(INC)
vpath %.ox.html $(DOC)

imports =AllTest.
includes=GetStarted.ox GetStartedData.ox HeteroData.ox BobsChoice.ox TommysChoice.ox 
exoxo = $(imports:.=.oxo)
exdoc = $(includes)  $(imports:.=.ox)

DDPexamples.oxo : $(includes) $(exoxo)

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	
.PHONY : document
document:
	$(OXDOC) -uplevel -include $(INC) $(exdoc)
