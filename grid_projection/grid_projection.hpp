#ifndef GRID_PROJ_H
#define GRID_PROJ_H

#include "../query_result/query_result.hpp"
#include "../hash_table/hash_table.hpp"

class Grid_Projection {
public:
	Grid_Projection(int L, int dimension, int w, int k, int delta,
			int curve_dimension, unsigned m);
	~Grid_Projection();

	void insert_curve(Curve *curve);
	void ANN(Curve *query_curve, unsigned threshhold, Query_Result& query_result,
		bool check_for_identical_grid_flag);
	void print_hash_tables();
	int get_w() {return w;}
	int get_k() {return k;}
	int get_dimension() {return dimension;}
	void snap_curve(Curve *curve, Point *t, Curve **snapped_curve, int delta);
	void get_snapped_point(Point *point, int delta, Point *t, Point **snapped_point);
	void zip_points(Curve *snapped_curve, Item **item);
	void fill_curve(Curve *curve, int pad_length);
	void convert_2d_curve_to_vector(Curve *curve, Point *t, int delta,
		Curve **snapped_curve, Item **item);
private:
	unsigned long long int Grid_Projection_distance(Curve *curve1, Curve *curve2);

	vector<Hash_Table*> hash_tables;
	vector<Point*> grids; //Each grid is identified from the random Point t
	int w;
	int k;
	int dimension;
	int L;
	int bits_of_each_hash;
	unsigned M;
	int delta;
	int curve_dimension; //2D in our case
	unsigned m;
	vector<unsigned> m_powers;
};

#endif
