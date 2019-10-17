CC = g++
CFLAGS = -std=c++11
OBJECTS = lsh_main.o lsh_implem.o hash_table_implem.o item_implem.o helping_functions.o
HEADERS = hash_table/hash_table.hpp item/item.hpp \
		helping_functions/helping_functions.hpp curve/curve.hpp grid_projection/grid_projection.hpp point/point.hpp
	
all: lsh_main

lsh_main: $(OBJECTS) 
	$(CC) $(OBJECTS) -o lsh_main 

lsh_main.o: lsh_main.cpp helping_functions/helping_functions.hpp
	$(CC) -c lsh_main.cpp -o lsh_main.o

lsh_implem.o: lsh/lsh_implem.cpp lsh/lsh.hpp 
	$(CC) -c lsh/lsh_implem.cpp -o lsh_implem.o 

hash_table_implem.o: hash_table/hash_table_implem.cpp hash_table/hash_table.hpp 
	$(CC) -c hash_table/hash_table_implem.cpp -o hash_table_implem.o 

item_implem.o: item/item_implem.cpp item/item.hpp
	$(CC) -c item/item_implem.cpp -o item_implem.o

helping_functions.o: helping_functions/helping_functions.cpp helping_functions/helping_functions.hpp 
	$(CC) -c helping_functions/helping_functions.cpp -o helping_functions.o

clean :
	rm *.o lsh_main