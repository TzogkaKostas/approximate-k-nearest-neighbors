#ifndef REL_TRAV_H
#define REL_TRAV_H

#include <list>
#include <vector>

#include "../query_result/query_result.hpp"
#include "../hash_table/hash_table.hpp"
#include "../Tuple/tuple.hpp"

class Relevant_Traversals {
public:
	Relevant_Traversals(int i, int j, int L, int K_matrix, int w, int k, int m, int M);
	Relevant_Traversals() {}
	~Relevant_Traversals();

	void insert(Curve *curve, int L, int w, int k, int bits_of_each_hash, int M,
		double **G_matrix, int K_matrix, int curve_dimensinion);
	void print_hash_tables();

	list<vector<Tuple*>*>& get_relevant_traversals() {return relevant_traversals;}
	Hash_Table** get_hash_tables() {return hash_tables;}
	vector<vector<unsigned>*>& get_m_powers_array() {return m_powers_array;}
	int get_num_of_traversals() {return relevant_traversals.size();}


private:
    list<vector<Tuple*>*> relevant_traversals;
	Hash_Table** hash_tables;
	vector<vector<unsigned>*> m_powers_array;
	int length_i;
	int length_j;
	int L;
};

#endif
