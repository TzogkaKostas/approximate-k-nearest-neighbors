#include "relevant_traversals_hypercube.hpp"
#include "../helping_functions//helping_functions.hpp"
#include "../curve_grid_hypercube/curve_grid_hypercube.hpp"
using namespace std;
//class Hash_Table_Hypercube;
Relevant_Traversals_hypercube::Relevant_Traversals_hypercube(int i, int j, int table_size_hypercube, int K_matrix, int w, int k,int m,int M){

    find_relevant_traversals(i + 1, j + 1, relevant_traversals);
    //cout <<"rela _size: ("<<i<<", "<<j<<"): "<<relevant_traversals.size()<<endl;

    int hash_table_dimension;
	hash_tables = new Hash_Table_Hypercube*[relevant_traversals.size()];
	int rel_indx = 0;
    for (vector<Tuple*>* rel_trav : relevant_traversals ) {
        hash_table_dimension = rel_trav->size()*K_matrix;
		
		hash_tables[rel_indx] = new Hash_Table_Hypercube(table_size_hypercube, hash_table_dimension, w, k);
		rel_indx++;

        vector<unsigned> *m_powers = new vector<unsigned>;
        for (size_t i = 0; i < hash_table_dimension; i++) {
            m_powers->push_back( pow_mod(m, i, M) );
        }
        m_powers_array.push_back(m_powers);
    }

    this->length_i = i;
    this->length_j = j;
    this->table_size_hypercube = table_size_hypercube;

    //cout <<"rela _size: ("<<i<<", "<<j<<"): "<<relevant_traversals.size()<<endl;

}

Relevant_Traversals_hypercube::~Relevant_Traversals_hypercube() {
	int i = 0;
	for (vector<Tuple*> *relevant_traversal : relevant_traversals) {
		for (Tuple *t : *relevant_traversal) {
			delete t;
		}
		delete relevant_traversal;
		delete hash_tables[i];
		delete m_powers_array[i];
		i++;
	}
	delete[] hash_tables;
}

void Relevant_Traversals_hypercube::insert(Curve *curve, int w,
    int k, int bits_of_each_hash, int M,
    double **G_matrix, int K_matrix, int curve_dimension) {

	Item *item = NULL;
	unsigned P_value;
    //cout<< "boom \n";
    //cout <<"\n rel size:\n"<<relevant_traversals.size()<<endl;
    //cout << "\n BOOM \n";
	int rel_indx = 0;
	for (vector<Tuple*> *relevant_traversal : relevant_traversals) {
		convert_2d_curve_to_vector_by_projection(*relevant_traversal, 0, G_matrix,
		curve, K_matrix, curve_dimension, &item);

		//cout <<"w: "<<w<<endl;
		//cout <<"k: "<<k<<endl;
		//cout <<"bits_of: "<<bits_of_each_hash<<endl;
		//cout <<"M: "<<M<<endl;
		//getchar();

		P_value = hash_tables[rel_indx]->p(*(item->get_coordinates()),hash_tables[rel_indx]->get_dimension(), table_size_hypercube, w,  k, bits_of_each_hash,  M, *m_powers_array[rel_indx]);
		//cout <<"p: "<<P_value<<endl;
		//getchar();
		hash_tables[rel_indx]->insert(curve, P_value);
		delete item;
		rel_indx++;

	}
}
