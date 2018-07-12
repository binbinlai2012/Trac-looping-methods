CXX = g++
CPPFLAGS = -O3
LDFLAGS = -O3

TARGET = ../bin/call_interaction_model ../bin/call_interaction ../bin/merge_sig_interaction ../bin/merge_tracloops ../bin/count_tag_on_region ../bin/get_petc_in_loop_from_bedpe ../bin/call_diff_loop
OBJECTS = operation.o pet.o partition.o call_interaction_model.o call_interaction.o merge_sig_interaction.o 

.PHONY : all
all : $(TARGET) 

../bin/call_interaction_model : operation.o pet.o call_interaction_model.o partition.o gene.o
	$(CXX) $(LDFLAGS) -o ../bin/call_interaction_model operation.o pet.o call_interaction_model.o partition.o gene.o

../bin/call_interaction : operation.o pet.o call_interaction.o
	$(CXX) $(LDFLAGS) -o ../bin/call_interaction operation.o pet.o call_interaction.o
	
../bin/merge_sig_interaction : operation.o pet.o merge_sig_interaction.o
	$(CXX) $(LDFLAGS) -o ../bin/merge_sig_interaction operation.o pet.o merge_sig_interaction.o
	
../bin/merge_tracloops : 
	$(CXX) $(LDFLAGS) -o ../bin/merge_tracloops merge_tracloops.cpp

../bin/count_tag_on_region :
	$(CXX) $(LDFLAGS) -o ../bin/count_tag_on_region count_tag_on_region.cpp

../bin/get_petc_in_loop_from_bedpe :
	$(CXX) $(LDFLAGS) -o ../bin/get_petc_in_loop_from_bedpe get_petc_in_loop_from_bedpe.cpp
	
../bin/call_diff_loop :
	$(CXX) $(LDFLAGS) -o ../bin/call_diff_loop call_diff_loop.cpp
	

operation.o :
pet.o : operation.h
gene.o : operation.h
partition.o : pet.h operation.h gene.h
call_interaction_model.o : pet.h partition.h
call_interaction.o : pet.h operation.h
merge_sig_interaction.o : pet.h operation.h




