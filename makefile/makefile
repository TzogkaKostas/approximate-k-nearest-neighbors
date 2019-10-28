CC           = g++
FLAGS        = -c
DEBUG        = -g3
OBJECTS =  curve_projection_hypercube.o point_implem.o item_implem.o helping_functions.o curve_implem.o curve_grid_hypercube.o relevant_traversals_hypercube.o
HEADERS =  ../hash_table/hash_table.hpp ../curve_projection_hypercube/curve_projection_hypercube.hpp ../relevant_traversals_hypercube/relevant_traversals_hypercube.hpp ../point/point.hpp  ../item/item.hpp ../query_result/query_result.hpp ../hypercube/hypercube.hpp ../Tuple/tuple.hpp ../helping_functions/helping_functions.hpp ../curve/curve.hpp ../curve_grid_hypercube/curve_grid_hypercube.hpp



#------------------------------------------------- #

all  : cube curve_grid_hypercube lsh curve_grid_lsh curve_projection_lsh

curve_implem.o : ../curve/curve_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../curve/curve_implem.cpp -o curve_implem.o

curve_grid_hypercube.o : ../curve_grid_hypercube/curve_grid_hypercube.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../curve_grid_hypercube/curve_grid_hypercube.cpp -o curve_grid_hypercube.o

curve_projection_lsh.o : ../curve_projection_lsh/curve_projection_lsh_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../curve_projection_lsh/curve_projection_lsh_implem.cpp -o curve_projection_lsh_implem.o

helping_functions.o : ../helping_functions/helping_functions.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../helping_functions/helping_functions.cpp -o helping_functions.o

hypercube.o : ../hypercube/hypercube.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../hypercube/hypercube.cpp -o hypercube.o

item_implem.o : ../item/item_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../item/item_implem.cpp -o item_implem.o

point_implem.o : ../point/point_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../point/point_implem.cpp -o point_implem.o

relevant_traversals.o : ../relevant_traversals/relevant_traversals_implem.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../relevant_traversals/relevant_traversals_implem.cpp -o relevant_traversals_implem.o

relevant_traversals_hypercube.o: ../relevant_traversals_hypercube/relevant_traversals_hypercube.cpp
	$(CC) $(DEBUG) $(FLAGS) ../relevant_traversals_hypercube/relevant_traversals_hypercube.cpp -o relevant_traversals_hypercube.o

curve_projection_hypercube.o : ../curve_projection_hypercube/curve_projection_hypercube.cpp $(HEADERS)
	$(CC) $(DEBUG) $(FLAGS) ../curve_projection_hypercube/curve_projection_hypercube.cpp -o curve_projection_hypercube.o

cube : ../hypercube_main/hypercube_main.cpp ../hypercube/hypercube.cpp $(OBJECTS)
	$(CC) $(DEBUG) $(OBJECTS) ../hypercube_main/hypercube_main.cpp ../hypercube/hypercube.cpp -o cube

curve_grid_hypercube : $(OBJECTS) ../curve_grid_hypercube_main/grid_hypercube_main.cpp
	$(CC) $(DEBUG) $(OBJECTS) ../curve_grid_hypercube_main/grid_hypercube_main.cpp -o curve_grid_hypercube

curve_projection_hypercube : $(OBJECTS) ../curve_projection_main/projection_hypercube_main.cpp
	$(CC) $(DEBUG) $(OBJECTS) ../curve_projection_main/projection_hypercube_main.cpp -o curve_projection_hypercube

lsh :  ../lsh_main/lsh_main.cpp ../hash_table_vector/hash_table_vector_implem.cpp ../lsh/lsh_implem.cpp $(OBJECTS)
	$(CC) $(DEBUG) $(OBJECTS) ../lsh_main/lsh_main.cpp ../hash_table_vector/hash_table_vector_implem.cpp ../lsh/lsh_implem.cpp -o lsh

curve_grid_lsh : ../grid_lsh_main/grid_lsh_main.cpp ../curve_grid_lsh/curve_grid_lsh_implem.cpp ../hash_table/hash_table_implem.cpp $(OBJECTS)
	$(CC) $(DEBUG) $(OBJECTS)	../grid_lsh_main/grid_lsh_main.cpp ../curve_grid_lsh/curve_grid_lsh_implem.cpp ../hash_table/hash_table_implem.cpp -o curve_grid_lsh

curve_projection_lsh : $(OBJECTS)  ../projection_lsh_main/projection_lsh_main.cpp ../curve_projection_lsh/curve_projection_lsh_implem.cpp ../relevant_traversals/relevant_traversals_implem.cpp ../hash_table/hash_table_implem.cpp
	$(CC) $(DEBUG) $(OBJECTS) ../projection_lsh_main/projection_lsh_main.cpp ../curve_projection_lsh/curve_projection_lsh_implem.cpp ../relevant_traversals/relevant_traversals_implem.cpp ../hash_table/hash_table_implem.cpp -o curve_projection_lsh

# Cleaning Files
clean:
	rm *.o cube curve_grid_hypercube curve_projection_hypercube lsh curve_grid_lsh curve_projection_lsh
