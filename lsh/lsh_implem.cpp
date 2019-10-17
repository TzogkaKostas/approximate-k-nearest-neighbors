#include <vector>
#include <list>
#include <cmath>
#include <iostream>
#include <string>
#include <limits>

using namespace std;

#include "../item/item.hpp"
#include "../hash_table/hash_table.hpp"
#include "../helping_functions/helping_functions.hpp"
#include "../query_result/query_result.hpp"
#include "lsh.hpp"

LSH::LSH(int L, int dimension, int w, int k, unsigned m) {
	for (size_t i = 0; i < L; i++) {
		Hash_Table *hash_table = new Hash_Table(dimension, w, k);
		hash_tables.push_back(hash_table);
	}
	this->L = L;
	this->dimension = dimension;
	this->w = w;    
	this->k = k;
	this->bits_of_each_hash = 32/k;
	if (k == 1) {
		this->M = numeric_limits<unsigned>::max();
	}
	else {
		this->M = pow(2, bits_of_each_hash);
	}
	this->m = m;
	for (size_t i = 0; i < dimension; i++) {
		m_powers.push_back( pow_mod(m, i, M) );
	}
}

LSH::~LSH() {
	for (size_t i = 0; i < L; i++) {
		delete hash_tables[i];
	}
}

void LSH::insert_item(Item *item) {
	unsigned g_value;
	for (size_t i = 0; i < L; i++) {
		g_value = g_hash_function(*(item->get_coordinates()), dimension,
			w, k, bits_of_each_hash, M, hash_tables[i]->get_s_array(), m_powers);
		hash_tables[i]->insert(item, g_value);
	}
}

void LSH::ANN(Item *query, unsigned threshhold, Query_Result& query_result) {
	unsigned searched_items;
	unsigned best_distance = numeric_limits<unsigned>::max();
	unsigned position, g_value;
	string best = "";
	unordered_multimap<unsigned, Item*> *map;
	pair <unordered_multimap<unsigned, Item*>::iterator, unordered_multimap<unsigned,Item*>::iterator> ret;
	unordered_multimap<unsigned, Item*>::iterator it;

	time_t time;
	time = clock();
	for (size_t i = 0; i < L; i++) {
		g_value = g_hash_function(*(query->get_coordinates()),
				dimension, w, k, bits_of_each_hash, M, hash_tables[i]->get_s_array(), m_powers);
		
		map = hash_tables[i]->get_map();
		ret = map->equal_range(g_value);
		searched_items = 0;
		for (it = ret.first; it != ret.second; ++it) {
			if (searched_items >= threshhold) {
				goto exit;
			}

			unsigned cur_distance = lsh_distance(query, it->second);
			if (cur_distance < best_distance) {
				best = it->second->get_name();
				best_distance = cur_distance;
			}
			searched_items++;
		}
	}
	exit:
	time = clock() - time;

	if (best != "") {
		query_result.set_best_distance(best_distance);
		query_result.set_time( ((double)time) / CLOCKS_PER_SEC );
		query_result.set_best_item(best);
	}
	else {
		query_result.set_best_distance(-1);
		query_result.set_time(-1);
		query_result.set_best_item("NULL");
	}
}

void LSH::print_hash_tables() {
	for (size_t i = 0; i < L; i++) {
		cout <<"Hash table: "<<i<<endl;
		hash_tables[i]->print();
	}
}

unsigned long long int LSH::lsh_distance(Item *item1, Item *item2) {
	return manhattan_distance(*(item1->get_coordinates()), *(item2->get_coordinates()));
}