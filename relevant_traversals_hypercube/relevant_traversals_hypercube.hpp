#ifndef REL_TRAV_HYPER_H
#define REL_TRAV_HYPER_H

#include <list>
#include <vector>

#include "../query_result/query_result.hpp"
#include "../Tuple/tuple.hpp"
#include "../curve_grid_hypercube/curve_grid_hypercube.hpp"


class Relevant_Traversals_hypercube {
public:
	Relevant_Traversals_hypercube(int i, int j, int table_size_hypercube, int K_matrix, int w, int k,int m,int M);
	~Relevant_Traversals_hypercube();

	void insert(Curve *curve, int w,
		int k, int bits_of_each_hash, int M,
		double **G_matrix, int K_matrix, int curve_dimensinion);
	void print_hash_tables();

	list<vector<Tuple*>*>& get_relevant_traversals() {return relevant_traversals;}
	Hash_Table_Hypercube** get_hash_tables() {return hash_tables;}
	vector<vector<unsigned>*>& get_m_powers_array() {return m_powers_array;}
	int get_num_of_traversals() {return relevant_traversals.size();}

private:
	list<vector<Tuple*>*> relevant_traversals;
	//vector<Hash_Table_Hypercube*>  hash_tables_R;
	Hash_Table_Hypercube** hash_tables;
	vector<vector<unsigned>*> m_powers_array;

	int table_size_hypercube;
	int length_i;
	int length_j;
};

#endif
