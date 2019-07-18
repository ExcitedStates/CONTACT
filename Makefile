ROOT = /Users/vdbedem/Downloads/mmdb-1.23.2.2
INCDIR = $(ROOT)/include/
LIBDIR = $(ROOT)/lib/
LIB = mmdb
CXX = g++ 
OBJECTS = mainfile.o CreateGraph_stats.o
ALL: GRAPH_STATS 
GRAPH_STATS: $(OBJECTS)
	$(CXX) -o Graph_STATS $(OBJECTS) $(addprefix -I, $(INCDIR)) $(addprefix -L, $(LIBDIR)) $(addprefix -l, $(LIB))
CreateGraph_stats.o: CreateGraph_stats.cpp
	$(CXX) -c CreateGraph_stats.cpp $(addprefix -I, $(INCDIR)) $(addprefix -L, $(LIBDIR)) $(addprefix -l, $(LIB))
mainfile.o: mainfile.cpp
	$(CXX) -c mainfile.cpp $(addprefix -I, $(INCDIR)) $(addprefix -L, $(LIBDIR)) $(addprefix -l, $(LIB))
.PHONY: clean
clean: 
	rm -f *.o *~
