CC = g++
CFLAGS = -std=c++11
OBJECTS = lsh_main.o lsh_implem.o hash_table_implem.o item_implem.o helping_functions.o \
		curve_implem.o point_implem.o
HEADERS = ../hash_table/hash_table.hpp ../item/item.hpp ../helping_functions/helping_functions.hpp \
		../lsh/lsh.hpp ../curve/curve.hpp ../point/point.hpp
	
all: lsh_main

lsh_main: $(OBJECTS) 
	$(CC) $(OBJECTS) -o lsh_main 

lsh_main.o: lsh_main.cpp $(HEADERS)
	$(CC) -c lsh_main.cpp -o lsh_main.o

lsh_implem.o: ../lsh/lsh_implem.cpp $(HEADERS) 
	$(CC) -c ../lsh/lsh_implem.cpp -o lsh_implem.o 

hash_table_implem.o: ../hash_table/hash_table_implem.cpp $(HEADERS)
	$(CC) -c ../hash_table/hash_table_implem.cpp -o hash_table_implem.o 

item_implem.o: ../item/item_implem.cpp $(HEADERS)
	$(CC) -c ../item/item_implem.cpp -o item_implem.o

helping_functions.o: ../helping_functions/helping_functions.cpp $(HEADERS) 
	$(CC) -c ../helping_functions/helping_functions.cpp -o helping_functions.o

curve_implem.o: ../curve/curve_implem.cpp $(HEADERS)
	$(CC) -c ../curve/curve_implem.cpp -o curve_implem.o

point_implem.o: ../point/point_implem.cpp $(HEADERS)
	$(CC) -c ../point/point_implem.cpp -o point_implem.o

clean :
	rm *.o lsh_main