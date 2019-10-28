#include "curve_projection_hypercube.hpp"

using namespace std;

Curve_Projection_hypercube::Curve_Projection_hypercube(int w, int k,
        int curve_dimension, unsigned m,unsigned M ,int K_matrix,int table_size_hypercube,int probes,int M_Table){

    this->K_matrix = K_matrix;
    this->curve_dimension=curve_dimension;
    this->w=w;
    this->k=k;
    this->bits_of_each_hash = 32/k;
    this->table_size= M_Table;
    if (k == 1) {
        this->M = numeric_limits<unsigned>::max() + 1;
    }
    else {
        this->M = pow(2, bits_of_each_hash);
    }
    this->m = m;
    this->table_size_hypercube=table_size_hypercube;

    //ALLOCATE M_Table*M_Table TABLE
    table = new Relevant_Traversals_hypercube**[table_size];
    for(int i = 0; i < table_size; ++i) {
        table[i] = new Relevant_Traversals_hypercube*[table_size];

        for(int j = 0; j < table_size; ++j) {
            table[i][j] = new Relevant_Traversals_hypercube(i, j,table_size_hypercube ,K_matrix, w, k, m, M);
        }

    //CREATE RANDOM G MATRIX ~N(0, 1)
    G_matrix = new double*[K_matrix];
    for (size_t i = 0; i < K_matrix; i++) {
        G_matrix[i] = new double[curve_dimension];
    }
    random_matrix(K_matrix, curve_dimension, G_matrix, 0, 1);
    }
}


Curve_Projection_hypercube::~Curve_Projection_hypercube() {
    for(int i = 0; i < table_size; ++i) {
	    for(int j = 0; j < table_size; ++j) {
			delete table[i][j];
		}
		delete[] table[i];
	}
	delete[] table;

	for (size_t i = 0; i < K_matrix; i++) {
		delete[] G_matrix[i];
	}

	delete[] G_matrix;
}


void Curve_Projection_hypercube::insert_curve(Curve *curve) {
	Curve *grid_curve = NULL;
	Item *item = NULL;
	unsigned g_value;

	int table_row = curve->get_length() - 1;
	//cout <<"\ntable_row\n: "<<table_row<<endl;
	for (size_t j = 0; j < table_size; j++) {
		//cout <<":-"<<table[table_row][j]->get_num_of_traversals() <<endl;
		table[table_row][j]->insert(curve, table_size_hypercube, w, k,
		bits_of_each_hash, M, G_matrix, K_matrix, curve_dimension);
	}
}

void Curve_Projection_hypercube::ANN(Curve *query_curve, unsigned probes, Query_Result& query_result) {
        unsigned searched_items;
    	double best_distance = numeric_limits<double>::max();
    	unsigned P_value;
    	Item *query_item;
    	string best = "";
    	unordered_multimap<unsigned, Curve*>::iterator it;
    	pair <unordered_multimap<unsigned, Curve*>::iterator,
    		unordered_multimap<unsigned,Curve*>::iterator> ret;

    	int table_column = query_curve->get_length() - 1;
    	int start_row = max(0, table_column - 2);
    	int end_row = min(table_size - 1, table_column + 2);
        int flag =-1;
    	time_t time;
    	time = clock();
        int bucket_value=0;
        for (size_t row = start_row; row < end_row; row++) {
            list<vector<Tuple*>*> relevant_traversals =
    			table[start_row][table_column]->get_relevant_traversals();

    		vector<Hash_Table_Hypercube*> hash_tables =
				table[start_row][table_column]->get_hash_tables();

    		vector<vector<unsigned>*> m_powers_array =
    			table[start_row][table_column]->get_m_powers_array();
            int h_i = 0;
            for (vector<Tuple*> *relevant_traversal : relevant_traversals) {
                convert_2d_curve_to_vector_by_projection(*relevant_traversal, 1, G_matrix,
          				query_curve, K_matrix, curve_dimension, &query_item);
                P_value = hash_tables[h_i]->p(*(query_item->get_coordinates()),hash_tables[h_i]->get_dimension(), table_size, w, k, bits_of_each_hash, M,  *m_powers_array[h_i]);
                int bucketes_checked=0;
                ret = hash_tables[h_i]->get_f_values_map()->equal_range(P_value);
        		searched_items = 0;
              for (it = ret.first; it != ret.second; ++it) {
        			if (searched_items >= M_f) {
        				flag =1;
        			}
                    if(flag == 1){
                        break;
                    }

        			unsigned cur_distance = Curve_Grid_distance(query_curve, it->second);
        			if (cur_distance < best_distance) {
        				best = it->second->get_name();
        				best_distance = cur_distance;
        			}
        			searched_items++;
        		}
                unsigned nbuckets=hash_tables[h_i]->get_f_values_map()->bucket_count();
                for (unsigned y=0; y<nbuckets; y++) {
                    if (searched_items >= M_f || bucketes_checked>=probes) {
        				break;
        			}
                    for (auto it = hash_tables[h_i]->get_f_values_map()->begin(y);it!= hash_tables[h_i]->get_f_values_map()->end(y);it++){
        	        	//std::cout << "[" << hash_table->get_f_values_map()->begin(i)->first << ":" << it->second->get_name() << "] ";
        				bucket_value=hash_tables[h_i]->get_f_values_map()->begin(y)->first;
        				break;
        	    	}
                    if (hammingDistance(P_value, bucket_value)==1){
        				bucketes_checked++;
        				ret = hash_tables[h_i]->get_f_values_map()->equal_range(bucket_value);
        				for (it = ret.first; it != ret.second; ++it) {
        					if (searched_items >= M_f) {
        						break;
        					}
        					unsigned cur_distance = Curve_Grid_distance((query_curve), (it->second));//apostasi querry apo ta alla pou iparxoun sto bucket
        					if (cur_distance < best_distance) {
        						best = it->second->get_name();
        						best_distance = cur_distance;
        					}
        					searched_items++;
        				}
        			}

                }
                h_i++;
            }
        }
}
double Curve_Projection_hypercube::Curve_Grid_distance(Curve *curve1, Curve *curve2) {
	return DTW(curve1->get_points(), curve2->get_points());
}
