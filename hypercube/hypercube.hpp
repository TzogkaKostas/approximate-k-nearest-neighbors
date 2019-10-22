#ifndef HYPERCUBE_H_
#define HYPERCUBE_H_
#include <unordered_map>
#include <math.h>
#include <limits>
#include <random>
#include <iostream>
#include "../query_result/query_result.hpp"
#include "../item/item.hpp"
#include "../helping_functions/helping_functions.hpp"

using namespace std;

typedef float Type;


class Hash_Table {
public:
	Hash_Table(int table_size, int dimension, int w, int k);
	~Hash_Table();
	void insert(Item*item, int dimension, int w, int k, int bits_of_each_has, unsigned M,vector<unsigned> m_powers);
	Item* find(Item *item, int dimension, int w, int k, int bits_of_each_has, unsigned M);
	unsigned p(vector<Type> x , int dimension, int table_size, int w, int k,
		int bits_of_each_hash, unsigned M, vector<unsigned>& m_powers);
	void print();
	int get_table_size() {return table_size;}
	vector<vector<float>*>& get_s_array(){return s_array;};
	unordered_multimap<unsigned, Item*>* get_f_values_map(){return &f_value;};

private:
	unordered_multimap <unsigned,Item*> f_value;
	vector< unordered_map<unsigned,int> *> g_value;
	vector < vector <float>* > s_array;
	int table_size;//k comnd line
};





class Hypercube{
public:
  	Hypercube(int hash_table_size, int dimension, int w, int k,unsigned m,unsigned M);
	void insert_item(Item *item);
	int get_w() {return w;}
	int get_k() {return k;}
	int get_dimension() {return dimension;}
	unsigned long long Hypercube_distance(Item *item1, Item *item2);
	void range_search(Item *query, unsigned threshhold, float radious,
		list<Item*>& range_items, Query_Result& query_result);
	void print_hash_tables();
	void ANN(Item *query, unsigned threshhold, Query_Result& query_result);
	Hash_Table* hash_table;
private:
	int table_size;//k comand line
  	int w;
  	int k;//for g,h
  	int dimension;
	int M_f;//M
  	unsigned long M;// for g,h
	unsigned m;
  	unsigned bits_of_each_hash;
	vector<unsigned> m_powers;
};
#endif
