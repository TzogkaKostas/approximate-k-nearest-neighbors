#include <vector>
#include <cmath>
#include <limits>

#include "../relevant_traversals/relevant_traversals.hpp"
#include "../helping_functions/helping_functions.hpp"


Relevant_Traversals::Relevant_Traversals(int i, int j, int L, int hash_table_dimension, int w, int k) {
	for (size_t i = 0; i < L; i++) {
		Hash_Table *hash_table = new Hash_Table(hash_table_dimension, w, k);
		hash_tables.push_back(hash_table);
	}
	this->length_i = i;
	this->length_j = j;
	find_relevant_traversals(i, j, relevant_traversals);
}

Relevant_Traversals::~Relevant_Traversals() {
	for (size_t i = 0; i < hash_tables.size(); i++) {
		delete hash_tables[i];
	}
}

void Relevant_Traversals::insert(Curve *curve, int hash_table_dimension, int w,
		int k, int bits_of_each_hash, int M, vector<unsigned>& m_powers,
		float **G_matrix, int K_matrix, int curve_dimension) {

	Item *item = NULL;
	unsigned g_value;

	int i = 0;
	for (vector<Tuple*> *relevant_traversal : relevant_traversals) {
		convert_2d_curve_to_vector_by_projection(*relevant_traversal, G_matrix,
			curve, K_matrix, curve_dimension, item);
		g_value = g_hash_function(*(item->get_coordinates()), hash_table_dimension,
			w, k, bits_of_each_hash, M, hash_tables[i]->get_s_array(), m_powers);
		hash_tables[i]->insert(curve, g_value);
		delete item;
		i++;
	}
}

void Relevant_Traversals::print_hash_tables() {
	for (size_t i = 0; i < hash_tables.size(); i++) {
		hash_tables[i]->print();
	}
}
