#ifndef HYPERCUBE_CURVED_GRIP_H_
#define HYPERCUBE_CURVED_GRIP_H_
#include <unordered_map>
#include <math.h>
#include <limits>
#include <random>
#include <iostream>
#include "../item/item.hpp"
#include "../curve/curve.hpp"
#include "../query_result/query_result.hpp"
#include "../point/point.hpp"
#include "../hash_table/hash_table.hpp"
#include "../helping_functions/helping_functions.hpp"
#include "../relevant_traversals/relevant_traversals.hpp"
#include "../Tuple/tuple.hpp"


typedef float Type;


class Hash_Table_Hypercube {
public:
	Hash_Table_Hypercube(int table_size, int dimension, int w, int k);
	~Hash_Table_Hypercube();
	void insert(Curve *curve, unsigned g_value);
	Curve* find(Curve *item, int dimension, int w, int k, int bits_of_each_has, unsigned M);
	unsigned p(vector<Type> x , int dimension, int table_size, int w, int k,
		int bits_of_each_hash, unsigned M, vector<unsigned>& m_powers);
	void print();
	int get_table_size() {return table_size;}
	int get_dimension(){return dimension;}
	//vector<vector<float>*>& get_s_array(){return s_array;};
	unordered_multimap<unsigned, Curve*>* get_f_values_map(){return &f_value;};

private:
	unordered_multimap <unsigned,Curve*> f_value;
	vector< unordered_map<unsigned,int> *> g_value;
	vector<vector < vector <float>* >* > s_array;
	int table_size;//k comnd line
	int dimension;
};


class Curve_Grid_hypercube {
 public:
    Curve_Grid_hypercube(int L, int hash_table_dimension, int w, int k, int delta,
			int curve_dimension, unsigned m,unsigned M ,int table_size,int probes);
    ~Curve_Grid_hypercube();
    void insert_curve(Curve *curve,list<Curve*> *grid_curves);
	int get_w() {return w;}
	int get_k() {return k;}
	int get_dimension() {return curve_dimension;}
	void print_hash_tables();
	void ANN(Curve *query_curve, unsigned threshhold, Query_Result& query_result,
			bool check_for_identical_grid_flag);
	unsigned long long int Curve_Grid_distance(Curve *curve1, Curve *curve2);
private:
    vector<Hash_Table_Hypercube*>  hash_tables;
    vector<Point*> grids; //Each grid is identified from the random Point t
	int table_size;//k comand line
    int L;
  	int w;
  	int k;//for g,h
	int hash_table_dimension;
	int M_f;//M
  	unsigned long M;// for g,h
	unsigned m;
  	unsigned bits_of_each_hash;
	vector<unsigned> m_powers;
    int delta;
	int curve_dimension; //2D in our case
};

#endif
