#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

#include "../relevant_traversals/relevant_traversals.hpp"
#include "../helping_functions/helping_functions.hpp"


Relevant_Traversals::Relevant_Traversals(int i, int j, int L, int K_matrix, int w, int k, int m, int M) {
	find_relevant_traversals(i + 1, j + 1, relevant_traversals);

	int hash_table_dimension;
	for (vector<Tuple*>* rel_trav : relevant_traversals ) {
		hash_table_dimension = rel_trav->size()*K_matrix;
		for (size_t i = 0; i < L; i++) {
			Hash_Table *hash_table = new Hash_Table(hash_table_dimension, w, k);
			hash_tables.push_back(hash_table);
		}
		vector<unsigned> *m_powers = new vector<unsigned>;
		for (size_t i = 0; i < hash_table_dimension; i++) {
			m_powers->push_back( pow_mod(m, i, M) );
		}
		m_powers_array.push_back(m_powers);
	}
	this->length_i = i;
	this->length_j = j;
	this->L = L;
}

Relevant_Traversals::~Relevant_Traversals() {
	int i = 0;
	for (vector<Tuple*> *relevant_traversal : relevant_traversals) {
		for (Tuple *t : *relevant_traversal) {
			delete t;
		}
		delete relevant_traversal;

		for (size_t j = 0; j < L; j++) {
			delete hash_tables[i*L + j];
		}
		delete m_powers_array[i];
		i++;
	}
}

void Relevant_Traversals::insert(Curve *curve, int hash_table_dimension, int L, int w,
		int k, int bits_of_each_hash, int M, float **G_matrix, int K_matrix,
		int curve_dimension) {

	Item *item = NULL;
	unsigned g_value;

	int rel_indx = 0;
	cout <<"kostas"<<endl;

	for (vector<Tuple*> *relevant_traversal : relevant_traversals) {
		for (size_t i = 0; i < L; i++) {
			convert_2d_curve_to_vector_by_projection(*relevant_traversal, G_matrix,
				curve, K_matrix, curve_dimension, item);
			cout <<"kostas"<<endl;
			
			g_value = g_hash_function(*(item->get_coordinates()), hash_table_dimension,
				w, k, bits_of_each_hash, M, hash_tables[rel_indx]->get_s_array(), *(m_powers_array[rel_indx]));
			hash_tables[rel_indx]->insert(curve, g_value);
			delete item;
			rel_indx++;
		}
	}
}

void Relevant_Traversals::print_hash_tables() {
	for (size_t i = 0; i < hash_tables.size(); i++) {
		hash_tables[i]->print();
	}
}
