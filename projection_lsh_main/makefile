CC = g++
CFLAGS = -std=c++11
OBJECTS = projection_lsh_main.o hash_table_implem.o item_implem.o \
		helping_functions.o curve_implem.o point_implem.o \
		relevant_traversals_implem.o curve_projection_lsh_implem.o

HEADERS = ../hash_table/hash_table.hpp ../item/item.hpp \
		../helping_functions/helping_functions.hpp ../curve/curve.hpp \
		../point/point.hpp ../query_result/query_result.hpp \
		../relevant_traversals/relevant_traversals.hpp ../Tuple/tuple.hpp \
		../curve_projection_lsh/curve_projection_lsh.hpp

all: projection_lsh_main

projection_lsh_main: $(OBJECTS)
	$(CC) $(OBJECTS) -o projection_lsh_main 

projection_lsh_main.o: projection_lsh_main.cpp $(HEADERS)
	$(CC) -c projection_lsh_main.cpp -o projection_lsh_main.o

hash_table_implem.o: ../hash_table/hash_table_implem.cpp $(HEADERS)
	$(CC) -c ../hash_table/hash_table_implem.cpp -o hash_table_implem.o 

item_implem.o: ../item/item_implem.cpp $(HEADERS)
	$(CC) -c ../item/item_implem.cpp -o item_implem.o

helping_functions.o: ../helping_functions/helping_functions.cpp $(HEADERS) 
	$(CC) -c ../helping_functions/helping_functions.cpp -o helping_functions.o

curve_implem.o: ../curve/curve_implem.cpp ../curve/curve.hpp ../point/point.hpp
	$(CC) -c ../curve/curve_implem.cpp -o curve_implem.o

point_implem.o: ../point/point_implem.cpp $(HEADERS)
	$(CC) -c ../point/point_implem.cpp -o point_implem.o

relevant_traversals_implem.o: ../relevant_traversals/relevant_traversals_implem.cpp $(HEADERS)
	$(CC) -c ../relevant_traversals/relevant_traversals_implem.cpp -o relevant_traversals_implem.o

curve_projection_lsh_implem.o: ../curve_projection_lsh/curve_projection_lsh_implem.cpp $(HEADERS)
	$(CC) -c ../curve_projection_lsh/curve_projection_lsh_implem.cpp -o curve_projection_lsh_implem.o

clean :
	rm *.o projection_lsh_main