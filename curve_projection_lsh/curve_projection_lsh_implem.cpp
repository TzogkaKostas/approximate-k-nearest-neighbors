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
#include "../Tuple/tuple.hpp"
#include "curve_projection_lsh.hpp"


Curve_Projection_LSH::Curve_Projection_LSH(int L, int w, int k,
	int curve_dimension, unsigned m, int M_table, int K_matrix) {

	this->K_matrix = K_matrix;
	this->table_size = M_table;
	this->L = L;
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

	//ALLOCATE M_Table*M_Table TABLE
	table = new Relevant_Traversals**[table_size];
	for(int i = 0; i < table_size; ++i) {
	    table[i] = new Relevant_Traversals*[table_size];
	    for(int j = 0; j < table_size; ++j) {
	        //table[i][j] = new Relevant_Traversals(i, j, L, K_matrix, w, k, m, M);
	        table[i][j] = NULL;
	    }
	}

	//CREATE RANDOM G MATRIX ~N(0, 1)
	G_matrix = new double*[K_matrix];
	for (size_t i = 0; i < K_matrix; i++) {
		G_matrix[i] = new double[curve_dimension];
	}
	random_matrix(K_matrix, curve_dimension, G_matrix, 0, 1);
}

Curve_Projection_LSH::~Curve_Projection_LSH() {
	for(int i = 0; i < table_size; ++i) {
	    for(int j = 0; j < table_size; ++j) {
			delete table[i][j];
		}
		delete[] table[i];
	}
	delete[] table;

	for (size_t i = 0; i < K_matrix; i++) {
		delete[] G_matrix[i];
	}

	delete[] G_matrix;
}

void Curve_Projection_LSH::insert_curve(Curve *curve) {
	Curve *grid_curve = NULL;
	Item *item = NULL;
	unsigned g_value;

	int table_row = curve->get_length() - 1;
	for (size_t j = 0; j < table_size; j++) {
		if (table[table_row][j] == NULL) {
	        table[table_row][j] = new Relevant_Traversals(table_row, j, L, K_matrix, w, k, m, M);
		}
		table[table_row][j]->insert(curve, L, w, k,
			bits_of_each_hash, M, G_matrix, K_matrix, curve_dimension);
	}
}

void Curve_Projection_LSH::ANN(Curve *query_curve, unsigned threshhold, Query_Result& query_result) {
	unsigned searched_items;
	double best_distance = numeric_limits<double>::max();
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
	for (size_t row = start_row; row < end_row; row++) {
		if (table[row][table_column] == NULL) {
			continue;
		}

		list<vector<Tuple*>*> relevant_traversals = 
			table[row][table_column]->get_relevant_traversals();

    	Hash_Table** hash_tables =
			table[row][table_column]->get_hash_tables();

		vector<vector<unsigned>*> m_powers_array = 
			table[row][table_column]->get_m_powers_array();


		int h_i = 0;
		for (vector<Tuple*> *relevant_traversal : relevant_traversals) {
			convert_2d_curve_to_vector_by_projection(*relevant_traversal, 1, G_matrix,
				query_curve, K_matrix, curve_dimension, &query_item);

			for (size_t j = 0; j < L; j++) {
				g_value = g_hash_function(*(query_item->get_coordinates()),
					hash_tables[h_i]->get_dimension(), w, k, bits_of_each_hash, M,
					hash_tables[h_i]->get_s_array(),(*m_powers_array[h_i/L]));

				ret = hash_tables[h_i]->get_map()->equal_range(g_value);
				searched_items = 0;
				for (it = ret.first; it != ret.second; ++it) {
					if (searched_items >= threshhold) {
						delete query_item;
						goto exit;
					}

					double cur_distance = Curve_Projection_LSH_distance(query_curve, it->second);
					if (cur_distance < best_distance) {
						best = it->second->get_name();
						best_distance = cur_distance;
					}
					searched_items++;
				}
				h_i++;
			}
			delete query_item;
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

double Curve_Projection_LSH::Curve_Projection_LSH_distance(Curve *curve1, Curve *curve2) {
	return DTW(curve1->get_points(), curve2->get_points());
}

void Curve_Projection_LSH::print_hash_tables() {
	int num_of_traversals = 0;
	for (size_t i = 0; i < table_size; i++) {
		for (size_t j = 0; j < table_size; j++) {
			cout <<table[i][j]->get_num_of_traversals()<<" ";
			num_of_traversals += table[i][j]->get_num_of_traversals();
		}
		cout <<endl;
	}
	cout <<"Total Number of Relevant traversals: "<<num_of_traversals<<endl;
}
