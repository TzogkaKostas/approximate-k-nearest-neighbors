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
#include "grid_projection.hpp"


Grid_Projection::Grid_Projection(int L, int hash_table_dimension, int w, int k, int delta,
		int curve_dimension, unsigned m) {

	for (size_t i = 0; i < L; i++) {
		Hash_Table *hash_table = new Hash_Table(hash_table_dimension, w, k);
		hash_tables.push_back(hash_table);
	}

	vector<float> *random_vector;
	Point *grid;

	//create a uniformly random vector t in [0,d)^d for each hash table
	for (size_t i = 0; i < L; i++) {
		random_vector = new vector<float>;
		random_float_vector(0, curve_dimension, *random_vector, curve_dimension);

		grid = new Point();
		for (float rand_coord: *random_vector) {
			grid->insert_coordinate(rand_coord);
		}
		grids.push_back(grid);
		delete random_vector;
	}

	this->L = L;
	this->hash_table_dimension = hash_table_dimension;
	this->curve_dimension = curve_dimension;
	this->w = w;    
	this->k = k;
	this->bits_of_each_hash = 32/k;
	this->delta = delta;
	if (k == 1) {
		this->M = numeric_limits<unsigned>::max();
	}
	else {
		this->M = pow(2, bits_of_each_hash);
	}
	this->m = m;
	for (size_t i = 0; i < hash_table_dimension; i++) {
		m_powers.push_back( pow_mod(m, i, M) );
	}
}

Grid_Projection::~Grid_Projection() {
	for (size_t i = 0; i < L; i++) {
		delete hash_tables[i];
		delete grids[i];
	}
}

void Grid_Projection::insert_curve(Curve *curve, list<Curve*> *grid_curves) {
	Curve *grid_curve = NULL;
	Item *item = NULL;
	unsigned g_value;

	for (size_t i = 0; i < L; i++) {
		convert_2d_curve_to_vector(curve, grids[i], delta, hash_table_dimension, curve_dimension,
			&grid_curve, &item);
		grid_curves->push_back(grid_curve);
		g_value = g_hash_function(*(item->get_coordinates()), hash_table_dimension,
			w, k, bits_of_each_hash, M, hash_tables[i]->get_s_array(), m_powers);
		hash_tables[i]->insert(grid_curve, g_value);
		delete item;
	}
}

void Grid_Projection::ANN(Curve *query_curve, unsigned threshhold, Query_Result& query_result,
		bool check_for_identical_grid_flag) {
	unsigned searched_items;
	unsigned best_distance = numeric_limits<unsigned>::max();
	unsigned position, g_value;
	Curve *query_grid_curve;
	Item *query_item;
	string best = "";
	pair <unordered_multimap<unsigned, Curve*>::iterator, unordered_multimap<unsigned,Curve*>::iterator> ret;
	unordered_multimap<unsigned, Curve*>::iterator it;

	time_t time;
	time = clock();
	for (size_t i = 0; i < L; i++) {
		//cout <<"curve: "<<endl;
		//query_curve->print();

		convert_2d_curve_to_vector(query_curve, grids[i], delta, hash_table_dimension,
		curve_dimension, &query_grid_curve, &query_item);
		//cout <<"grid curve: "<<endl;
		//query_grid_curve->print();
		//cout <<"item: "<<endl;
		//query_item->print();
		g_value = g_hash_function(*(query_item->get_coordinates()),
				hash_table_dimension, w, k, bits_of_each_hash, M, hash_tables[i]->get_s_array(), m_powers);
		
		//cout <<"g_value: "<< g_value<<endl;

		ret = hash_tables[i]->get_map()->equal_range(g_value);
		searched_items = 0;
		for (it = ret.first; it != ret.second; ++it) {
			if (searched_items >= threshhold) {
				goto exit;
			}
			if (check_for_identical_grid_flag == true) {
				if (it->second->get_corresponding_curve()->identical(query_curve) == false) {
					continue;
				}
			}

			unsigned cur_distance = Grid_Projection_distance(query_curve, it->second);
			if (cur_distance < best_distance) {
				best = it->second->get_name();
				best_distance = cur_distance;
			}
			searched_items++;
		}
		delete query_grid_curve;
		delete query_item;
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

unsigned long long int Grid_Projection::Grid_Projection_distance(Curve *curve1, Curve *curve2) {
	return DTW(curve1->get_points(), curve2->get_points());
}

void Grid_Projection::print_hash_tables() {
	for (size_t i = 0; i < L; i++) {
		cout <<"Hash table: "<<i<<endl;
		hash_tables[i]->print();
	}
}