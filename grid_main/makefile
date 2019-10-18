CC = g++
CFLAGS = -std=c++11
OBJECTS = grid_main.o hash_table_implem.o item_implem.o \
		helping_functions.o curve_implem.o grid_projection_implem.o point_implem.o

HEADERS = ../hash_table/hash_table.hpp ../item/item.hpp \
		../helping_functions/helping_functions.hpp ../curve/curve.hpp ../grid_projection/grid_projection.hpp \
		../point/point.hpp ../query_result/query_result.hpp

all: grid_main

grid_main: $(OBJECTS)
	$(CC) $(OBJECTS) -o grid_main 

grid_main.o: grid_main.cpp $(HEADERS)
	$(CC) -c grid_main.cpp -o grid_main.o

hash_table_implem.o: ../hash_table/hash_table_implem.cpp $(HEADERS)
	$(CC) -c ../hash_table/hash_table_implem.cpp -o hash_table_implem.o 

item_implem.o: ../item/item_implem.cpp $(HEADERS)
	$(CC) -c ../item/item_implem.cpp -o item_implem.o

helping_functions.o: ../helping_functions/helping_functions.cpp $(HEADERS) 
	$(CC) -c ../helping_functions/helping_functions.cpp -o helping_functions.o

curve_implem.o: ../curve/curve_implem.cpp ../curve/curve.hpp ../point/point.hpp
	$(CC) -c ../curve/curve_implem.cpp -o curve_implem.o

grid_projection_implem.o: ../grid_projection/grid_projection_implem.cpp $(HEADERS)
	$(CC) -c ../grid_projection/grid_projection_implem.cpp -o grid_projection_implem.o

point_implem.o: ../point/point_implem.cpp $(HEADERS)
	$(CC) -c ../point/point_implem.cpp -o point_implem.o

clean :
	rm *.o grid_main