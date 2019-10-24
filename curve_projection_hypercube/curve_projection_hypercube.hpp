#ifndef HYPERCUBE_CURVED_PROJECTION_H_
#define HYPERCUBE_CURVED_PROJECTION_H_
#include <unordered_map>
#include <math.h>
#include <limits>
#include <random>
#include <iostream>
#include "../item/item.hpp"
#include "../curve/curve.hpp"
#include "../query_result/query_result.hpp"
#include "../point/point.hpp"
#include "../helping_functions/helping_functions.hpp"
#include "../relevant_traversals_hypercube/relevant_traversals_hypercube.hpp"
#include "../Tuple/tuple.hpp"

class Curve_Projection_hypercube {
 public:
    Curve_Projection_hypercube(int hash_table_dimension, int w, int k, int delta,
			int curve_dimension, unsigned m,unsigned M ,int table_size,int probes);
    ~Curve_Projection_hypercube();
    void insert_curve(Curve *curve,list<Curve*> *grid_curves);
	int get_w() {return w;}
	int get_k() {return k;}
	int get_dimension() {return curve_dimension;}
	void print_hash_tables();
	void ANN(Curve *query_curve, unsigned threshhold, Query_Result& query_result,
			bool check_for_identical_grid_flag);
	unsigned long long int Curve_Grid_distance(Curve *curve1, Curve *curve2);
private:

    vector<Point*> grids; //Each grid is identified from the random Point t
	int table_size;//k comand line
  	int w;
  	int k;//for g,h
	int hash_table_dimension;
	int M_f;//M
  	unsigned long M;// for g,h
	unsigned m;
  	unsigned bits_of_each_hash;
	vector<unsigned> m_powers;
	int curve_dimension; //2D in our case
    vector<vector<Relevant_Traversals*>> table;
    float **G_matrix;
	int K_matrix;
};
#endif
