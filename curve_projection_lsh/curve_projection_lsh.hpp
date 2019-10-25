#ifndef PROJ_H
#define PROJ_H

#include "../query_result/query_result.hpp"
#include "../hash_table/hash_table.hpp"
#include "../relevant_traversals/relevant_traversals.hpp"

class Curve_Projection_LSH {
public:
	Curve_Projection_LSH(int L, int w, int k, int curve_dimension, unsigned m,
		int max_curve_length, int K_matrix);
	~Curve_Projection_LSH();

	void insert_curve(Curve *curve);
	void ANN(Curve *query_curve, unsigned threshhold, Query_Result& query_result);
	void print_hash_tables();
	int get_w() {return w;}
	int get_k() {return k;}
private:
	double Curve_Projection_LSH_distance(Curve *curve1, Curve *curve2);
	int w;
	int k;
	int L;
	int bits_of_each_hash;
	unsigned long M;
	int curve_dimension; //2D in our case
	unsigned m;
	int table_size;
	Relevant_Traversals ***table;
	float **G_matrix;
	int K_matrix;
};

#endif
