include ../../source/Makefile.inc
DOC := ../../docs/examples/CFMPI
DOC := ..\..\docs\examples\CFMPI
INC := ../../include

vpath %.h .:$(INC)
vpath %.oxo .:$(INC)
vpath %.ox.html $(DOC)

imports = baseMPI. ClientServerTestA. PeerCommunicationTestA. 
exoxo = $(imports:.=.oxo)
exdoc =  $(imports:.=.ox)

MPItests.oxo : $(exoxo)

%.oxo : %.ox %.h
	$(OX) $(OXFLAGS) -i$(INC) -c $<
	
.PHONY : clean
clean:
	$(info $(CURDIR))
	$(ERASE) /Q *.oxo

.PHONY : document
document:
	$(OXDOC) -uplevel -include $(INC) $(exdoc)
	$(ERASE) /Q $(MDOC)\*.css
