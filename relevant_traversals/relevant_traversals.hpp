#ifndef REL_TRAV_H
#define REL_TRAV_H

#include <list>
#include <vector>

#include "../query_result/query_result.hpp"
#include "../hash_table/hash_table.hpp"
#include "../Tuple/tuple.hpp"

class Relevant_Traversals {
public:
	Relevant_Traversals(int i, int j, int L, int hash_table_dimension, int w, int k);
	~Relevant_Traversals();

	void insert(Curve *curve, int hash_table_dimension, int w,
		int k, int bits_of_each_hash, int M, vector<unsigned>& m_powers,
		float **G_matrix, int K_matrix, int curve_dimensinion);
	void print_hash_tables();

	list<vector<Tuple*>*>& get_relevant_traversals() {return relevant_traversals;}
	vector<Hash_Table*>& get_hash_tables() {return hash_tables;}


private:
    list<vector<Tuple*>*> relevant_traversals;
    vector<Hash_Table*> hash_tables;

	int length_i;
	int length_j;
};

#endif
