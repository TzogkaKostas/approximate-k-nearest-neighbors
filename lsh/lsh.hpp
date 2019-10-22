#ifndef LSH_H_
#define LSH_H_

#include "../query_result/query_result.hpp"
#include "../hash_table/hash_table.hpp"
#include <list>
using namespace std;
class LSH {
public:
	LSH(int L, int dimension, int w, int k, unsigned m);
	~LSH();

	void insert_item(Item *item);
	void ANN(Item *query, unsigned threshhold, Query_Result& query_result);
	void range_search(Item *query, unsigned threshhold, float radious,
		list<Item*>& range_items, Query_Result& query_result);
	void print_hash_tables();
	int get_w() {return w;}
	int get_k() {return k;}
	int get_dimension() {return dimension;}
private:
	unsigned long long int lsh_distance(Item *item1, Item *item2);
	vector<Hash_Table*> hash_tables;
	int w;
	int k;
	int dimension;
	int L;
	int bits_of_each_hash;
	unsigned long M;
	unsigned m;
	vector<unsigned> m_powers;
};
#endif
