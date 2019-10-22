#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <cmath>

using namespace std;

#include "../item/item.hpp"
#include "../curve/curve.hpp"
#include "../query_result/query_result.hpp"
#include "../point/point.hpp"
#include "../hash_table/hash_table.hpp"
#include "../helping_functions/helping_functions.hpp"
#include "../relevant_traversals/relevant_traversals.hpp"
#include "../tuple/tuple.hpp"
#include "curve_projection_lsh.hpp"


Curve_Projection_LSH::Curve_Projection_LSH(int L, int hash_table_dimension,
	int w, int k, int curve_dimension, unsigned m, int max_curve_length, int K_matrix) {

	//ALLOCATE M*M TABLE
	Relevant_Traversals *relevant_traversals;
	for (size_t i = 0; i < max_curve_length; i++) {
		for (size_t j = 0; j < max_curve_length; j++) {
			relevant_traversals = new Relevant_Traversals(i, j, L, hash_table_dimension, w, k);
			table[i].push_back(relevant_traversals);
		}
	}

	//CREATE RANDOM G MATRIX ~N(0, 1)
	random_matrix(K_matrix, curve_dimension, G_matrix, 0, 1);
	
	this->K_matrix = K_matrix;
	this->table_size = max_curve_length;
	this->L = L;
	this->hash_table_dimension = hash_table_dimension;
	this->curve_dimension = curve_dimension;
	this->w = w;    
	this->k = k;
	this->bits_of_each_hash = 32/k;
	if (k == 1) {
		this->M = numeric_limits<unsigned>::max() + 1;
	}
	else {
		this->M = pow(2, bits_of_each_hash);
	}
	this->m = m;
	for (size_t i = 0; i < hash_table_dimension; i++) {
		m_powers.push_back( pow_mod(m, i, M) );
	}
}

Curve_Projection_LSH::~Curve_Projection_LSH() {
	for (size_t i = 0; i < L; i++) {
		for (size_t j = 0; j < L; j++) {
			delete table[i][j];
		}
	}
}

void Curve_Projection_LSH::insert_curve(Curve *curve) {
	Curve *grid_curve = NULL;
	Item *item = NULL;
	unsigned g_value;

	int table_row = curve->get_length();
	for (size_t j = 0; j < table_size; j++) {
		table[table_row][j]->insert(curve, hash_table_dimension, w, k,
			bits_of_each_hash, M, m_powers, G_matrix, table_size, curve_dimension);
	}
}

void Curve_Projection_LSH::ANN(Curve *query_curve, unsigned threshhold, Query_Result& query_result) {
	unsigned searched_items;
	unsigned best_distance = numeric_limits<unsigned>::max();
	unsigned g_value;
	Item *query_item;
	string best = "";
	unordered_multimap<unsigned, Curve*>::iterator it;
	pair <unordered_multimap<unsigned, Curve*>::iterator,
		unordered_multimap<unsigned,Curve*>::iterator> ret;

	int table_column = query_curve->get_length() - 1;
	int start_row = max(0, table_column - 2);
	int end_row = min(table_size - 1, table_column + 2);

	time_t time;
	time = clock();
	for (size_t i = start_row; i < end_row; i++) {
		list<vector<Tuple*>*> relevant_traversals = 
			table[start_row][table_column]->get_relevant_traversals();

    	vector<Hash_Table*> hash_tables =
			table[start_row][table_column]->get_hash_tables();

		int h_i = 0;
		for (vector<Tuple*> *relevant_traversal : relevant_traversals) {
			convert_2d_curve_to_vector_by_projection(*relevant_traversal, G_matrix,
				query_curve, K_matrix, curve_dimension, query_item);
			g_value = g_hash_function(*(query_item->get_coordinates()), hash_table_dimension,
				w, k, bits_of_each_hash, M, hash_tables[i]->get_s_array(), m_powers);

			
			ret = hash_tables[i]->get_map()->equal_range(g_value);
			searched_items = 0;
			for (it = ret.first; it != ret.second; ++it) {
				if (searched_items >= threshhold) {
					goto exit;
				}

				unsigned cur_distance = Curve_Projection_LSH_distance(query_curve, it->second);
				if (cur_distance < best_distance) {
					best = it->second->get_name();
					best_distance = cur_distance;
				}
				searched_items++;
			}
			delete query_item;
			h_i++;
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

unsigned long long int Curve_Projection_LSH::Curve_Projection_LSH_distance(Curve *curve1, Curve *curve2) {
	return DTW(curve1->get_points(), curve2->get_points());
}

void Curve_Projection_LSH::print_hash_tables() {
	for (size_t i = 0; i < table_size; i++) {
		for (size_t j = 0; j < table_size; j++) {
			table[i][j]->print_hash_tables();\
		}
	}
}