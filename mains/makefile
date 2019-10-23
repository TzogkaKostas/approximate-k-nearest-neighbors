CC           = g++
FLAGS        = -c
DEBUG        = -g3
OBJECTS =Tuple.o relevant_traversals.o query_result.o point.o lsh.o item.o hypercube.o helping_functions.o curve.o curve_grid_hypercube.o curve_grid_lsh.o curve_projection_lsh.o hash_table.o
HEADERS =  relevant_traversals/relevant_traversals.hpp ../point/point.hpp ../lsh/lsh.hpp ../item/item.hpp ../hypercube/hypercube.hpp ../hash_table/hash_table_implem.hpp ../helping_functions/helping_functions.hpp ../curve_projection_lsh/curve_projection_lsh.hpp ../curve_grid_lsh/curve_grid_lsh.hpp ../curve/curve_implem.hpp ../curve_grid_hypercube/curve_grid_hypercube.hpp
#------------------------------------------------- #

all  : cube

curve.o : ../curve/curve_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../curve/curve_implem.cpp -o

curve_grid_hypercube.o : ../curve_grid_hypercube/curve_grid_hypercube.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../curve_grid_hypercube/curve_grid_hypercube.cpp -o curve_grid_hypercube.o

curve_grid_lsh.o : ../curve_grid_lsh/curve_grid_lsh_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS)	../curve_grid_lsh/curve_grid_lsh_implem.cpp -o curve_grid_lsh_implem.o

curve_projection_lsh.o : ../curve_projection_lsh/curve_projection_lsh_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../curve_projection_lsh/curve_projection_lsh_implem.cpp -o curve_projection_lsh_implem.o

hash_table.o : ../hash_table/hash_table_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../hash_table/hash_table_implem.cpp -o hash_table_implem.o

helping_functions.o : ../helping_functions/helping_functions.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../helping_functions/helping_functions.cpp -o helping_functions.o

hypercube.o : ../hypercube/hypercube.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../hypercube/hypercube.cpp -o hypercube

item.o : ../item/item_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../item/item_implem.cpp -o item_implem.o

lsh.o : ../lsh/lsh_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../lsh/lsh_implem.cpp -o lsh_implem.o

point.o : ../point/point_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../point/point_implem.cpp -o point_implem.o

query_result.o : query_result/query_result.hpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) query_result/query_result.hpp -o query_result.o

relevant_traversals.o : relevant_traversals/relevant_traversals_implem.cpp $(HEADERS)
 	$(CC) $(DEBUG) $(FLAGS) relevant_traversals/relevant_traversals_implem.cpp -o relevant_traversals_implem.o

Tuple.o : Tuple/tuple.hpp $(HEADERS)
 	$(CC) $(DEBUG) $(FLAGS) Tuple/tuple.hpp -o Tuple/tuple.o

cube_main.o : hypercube_main.o $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) hypercube_main.cpp -o hypercube_main.o

cube : $(OBJECTS)
	$(CC) $(DEBUG) $(OBJECTS) -o cube

# Cleaning Files
clean:
	rm *.o
