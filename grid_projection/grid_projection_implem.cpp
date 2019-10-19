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


Grid_Projection::Grid_Projection(int L, int dimension, int w, int k, int delta,
		int curve_dimension, unsigned m) {

	for (size_t i = 0; i < L; i++) {
		Hash_Table *hash_table = new Hash_Table(dimension, w, k);
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
	this->dimension = dimension;
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
	for (size_t i = 0; i < dimension; i++) {
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
		convert_2d_curve_to_vector(curve, grids[i], delta, &grid_curve, &item);
		grid_curves->push_back(grid_curve);
		g_value = g_hash_function(*(item->get_coordinates()), dimension,
			w, k, bits_of_each_hash, M, hash_tables[i]->get_s_array(), m_powers);
		hash_tables[i]->insert(grid_curve, g_value);
		delete item;
	}
}

void Grid_Projection::convert_2d_curve_to_vector(Curve *curve, Point *t, int delta,
		Curve **grid_curve, Item **item) {
	
	snap_curve(curve, t, grid_curve, delta);
	fill_curve(*grid_curve, dimension/curve_dimension - (*grid_curve)->get_length()); //padding
	(*grid_curve)->set_corresponding_curve(curve);
	zip_points(*grid_curve, item);
}

void Grid_Projection::snap_curve(Curve *curve, Point *t, Curve **grid_curve, int delta) {

	Point *snapped_point = NULL;

	vector<Point*> *snapped_points = new vector<Point*>;
	for(Point *point : curve->get_points() ) {
		get_snapped_point(point, delta, t, &snapped_point);
		//consecutive duplicate
		if (snapped_points->size() > 0) {
			if (snapped_points->back()->equals(snapped_point) == true) {
				delete snapped_point;
				continue;
			}
		}
		snapped_points->push_back(snapped_point);
	}
	*grid_curve = new Curve(snapped_points);

}

void Grid_Projection::get_snapped_point(Point *point, int delta, Point *t, Point **snapped_point) {
	/*
			one cell, shifted by t, of the big Grid
		up_left -------------- up_right
				|\ d1	 d2 /|
				| \		   / |
				|	 point   |
				|  / d3	  \d4|
		down_left ------------- down_left
	*/
	//cout <<"t_x: "<<t->get_x()<<endl;
	//cout <<"t_y: "<<t->get_y()<<endl;
	//cout <<"point_x: "<<point->get_x()/delta<<endl;
	//cout <<"point_y: "<<point->get_y()/delta<<endl;
	//cout <<"delta: "<<delta<<endl;


	float up_left_x = floor( point->get_x()/delta )*delta + t->get_x();
	float up_left_y = ceil( point->get_y()/delta )*delta + t->get_y();
	//cout <<"up_left_x: "<<up_left_x<<endl;
	//cout <<"up_left_y: "<<up_left_y<<endl;

	float up_right_x = ceil( point->get_x()/delta)*delta + t->get_x();
	float up_right_y = ceil( point->get_y()/delta)*delta + t->get_y();
	//cout <<"up_right_x: "<<up_right_x<<endl;
	//cout <<"up_right_y: "<<up_right_y<<endl;

	float down_left_x = floor( point->get_x()/delta)*delta + t->get_x();
	float down_left_y = floor( point->get_y()/delta)*delta + t->get_y();

	//cout <<"down_left_x: "<<down_left_x<<endl;
	//cout <<"down_left_y: "<<down_left_y<<endl;

	float down_right_x = ceil( point->get_x()/delta)*delta + t->get_x();
	float down_right_y = floor( point->get_y()/delta)*delta + t->get_y();

	//cout <<"down_right_x: "<<down_right_x<<endl;
	//cout <<"down_right_y: "<<down_right_y<<endl;

	//maxhattan distances 
	float d1 = abs(up_left_x - point->get_x()) + abs(up_left_y - point->get_y());
	float d2 = abs(up_right_x - point->get_x()) + abs(up_right_y - point->get_y());
	float d3 = abs(down_left_x - point->get_x()) + abs(down_left_y - point->get_y());
	float d4 = abs(down_right_x - point->get_x()) + abs(down_right_y - point->get_y());

	float min_d = d1;
	float best_x = up_left_x;
	float best_y = up_left_y;

	if (d2 < min_d) {
		min_d = d2;
		best_x = up_right_x;
		best_y = up_right_y;
	}

	if (d3 < min_d) {
		min_d = d3;
		best_x = down_left_x;
		best_y = down_left_y;
	}
	
	if (d4 < min_d) {
		min_d = d4;
		best_x = down_right_x;
		best_y = down_right_y;
	}
	*snapped_point = new Point(best_x, best_y);
}

void Grid_Projection::fill_curve(Curve *curve, int pad_length) {
	for (size_t i = 0; i < pad_length; i++) {
		Point *point = new Point(0, 0);
		curve->insert_point(point);
	}
}

void Grid_Projection::zip_points(Curve *grid_curve, Item **item) {
	vector<Type> *coordinates = new vector<Type>;
	for(Point *point : grid_curve->get_points() ) {
		coordinates->push_back(point->get_x());
		coordinates->push_back(point->get_y());
	}
	*item = new Item(coordinates);
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

		convert_2d_curve_to_vector(query_curve, grids[i], delta, &query_grid_curve, &query_item);
		//cout <<"grid curve: "<<endl;
		//query_grid_curve->print();
		//cout <<"item: "<<endl;
		//query_item->print();
		g_value = g_hash_function(*(query_item->get_coordinates()),
				dimension, w, k, bits_of_each_hash, M, hash_tables[i]->get_s_array(), m_powers);
		
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