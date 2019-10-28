#ifndef GRID_PROJ_H
#define GRID_PROJ_H

#include "../query_result/query_result.hpp"
#include "../hash_table/hash_table.hpp"
#include <list>

class Curve_Grid_LSH {
public:
	Curve_Grid_LSH(int L, int hash_table_dimension, int w, int k, float delta,
			int curve_dimension, unsigned m, double max_coord);
	~Curve_Grid_LSH();

	void insert_curve(Curve *curve, list<Curve*> *grid_curves);
	void ANN(Curve *query_curve, unsigned threshhold, Query_Result& query_result,
		bool check_for_identical_grid_flag);
	void print_hash_tables();
	void print_hash_tables_names();
	int get_w() {return w;}
	int get_k() {return k;}
	int get_dimension() {return hash_table_dimension;}

private:
	double Curve_Grid_LSH_distance(Curve *curve1, Curve *curve2);

	vector<Hash_Table*> hash_tables;
	vector<Point*> grids; //Each grid is identified from the random Point t
	int w;
	int k;
	int hash_table_dimension;
	int L;
	int bits_of_each_hash;
	unsigned long M;
	float delta;
	int curve_dimension; //2D in our case
	unsigned m;
	double max_coord;
	vector<unsigned> m_powers;
};

#endif
