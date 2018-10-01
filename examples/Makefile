include ../source/Makefile.inc
DOC := ../docs/examples
INC := ../include

vpath %.ox .:replications:DDP:FiveO:CFMPI:misc
vpath %.h .:$(INC):replications:DDP:FiveO:CFMPI:misc
vpath %.oxo .:$(INC):replications:DDP:FiveO:CFMPI:misc
vpath %.ox.html $(DOC)

imports = DDP/examples. FiveO/examples. replications. MPItests. misc/examples.
exoxo = $(imports:.=.oxo)
exdoc = main.ox 

#DDP/DDPexamples.ox FiveO/FiveOexamples.ox replications/replications.ox CFMPI/MPItests.ox misc/MiscExamples.ox

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
	${MAKE} -C DDP document
	${MAKE} -C FiveO document
	${MAKE} -C misc document
	${MAKE} -C CFMPI document
	${MAKE} -C replications document
