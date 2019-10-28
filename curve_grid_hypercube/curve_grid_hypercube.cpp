#include "curve_grid_hypercube.hpp"

using namespace std;

Curve_Grid_hypercube::Curve_Grid_hypercube(int L, int hash_table_dimension, int w, int k, double delta,
        int curve_dimension, unsigned m,unsigned M ,int table_size,int probes, double max_coord){
    //printf("\n %d %d %d %d %d %d %ld %d \n", L,hash_table_size,  curve_dimension, w, k, delta,m,M);
    for (int i = 0; i < L; i++) {
        Hash_Table_Hypercube *hash_table = new Hash_Table_Hypercube(table_size,hash_table_dimension ,w, k);
        hash_tables.push_back(hash_table);
    }
    vector<float> *random_vector;
    Point *grid;
    //create a uniformly random vector t in [0,d)^d for each hash table
    for (int i = 0; i < L; i++) {
        random_vector = new vector<float>;
        random_float_vector(0, delta, *random_vector, curve_dimension);


        grid = new Point();
        for (float rand_coord: *random_vector) {
            grid->insert_coordinate(rand_coord);
            //cout << "boom \n";
        }
        grids.push_back(grid);
        //grid->print_coordinates();
        //cout <<endl;
        delete random_vector;
    }

  this->max_coord=max_coord;
  this->L = L;
	this->table_size = table_size;
  this->hash_table_dimension = hash_table_dimension;
	this->curve_dimension = curve_dimension;
	this->w = w;
	this->k = k;
    this->M_f=M;
	this->bits_of_each_hash = 32/k;
	this->delta = delta;
	if (k == 1) {
		this->M = numeric_limits<unsigned>::max() ;
	}
	else {
		this->M = pow(2, bits_of_each_hash);
	}
	this->m = m;
    //cout << ":table_size "<<endl;
	for (size_t i = 0; i < table_size; i++) {
		//cout <<"m: "<<m<<endl;
		//cout <<"i: "<<i<<endl;
		//cout <<"M: "<<M<<endl;
		m_powers.push_back( pow_mod(m, i, M) );
		//cout << pow(m, i)<<endl;
		//cout << m_powers[i]<<endl;
		//getchar();
        //cout << "boom \n";
        //cout << "boom \n";
	}

}

Curve_Grid_hypercube::~Curve_Grid_hypercube(){
    for (size_t i = 0; i < L; i++) {
		delete hash_tables[i];
		delete grids[i];
    }
}

void Curve_Grid_hypercube::insert_curve(Curve *curve, list<Curve*> *grid_curves) {
	Curve *grid_curve = NULL;
	Item *item = NULL;
	unsigned P_value;

    for (size_t i = 0; i < L; i++) {
        //curve->print();
        //grids[i]->print_coordinates();
        //cout <<endl <<delta <<endl<< table_size<<endl<<curve_dimension<<endl;
        //cout << endl;

        convert_2d_curve_to_vector(curve, grids[i], delta, hash_table_dimension, curve_dimension,
			&grid_curve, &item,max_coord);
        //cout << "\n input"<<endl;
		grid_curves->push_back(grid_curve);
		P_value = hash_tables[i]->p(*item->get_coordinates(), hash_table_dimension, table_size, w, k,bits_of_each_hash,  M, m_powers);
		hash_tables[i]->insert(grid_curve, P_value);
		//cout <<"P_value: "<<P_value<<endl;
		delete item;
	}
}

void Curve_Grid_hypercube::ANN(Curve *query_curve, unsigned probes, Query_Result& query_result,bool check_for_identical_grid_flag) {
    unsigned searched_items;
    double best_distance = numeric_limits<double>::max();
    unsigned P_value;
    Curve *query_grid_curve;
    Item *query_item;
    string best = "";
    pair <unordered_multimap<unsigned, Curve*>::iterator, unordered_multimap<unsigned,Curve*>::iterator> ret;
    unordered_multimap<unsigned, Curve*>::iterator it;
    int flag =-1;
    time_t time;
    unsigned bucket_value;
    time = clock();
    for (size_t i = 0; i < L; i++) {
        //cout <<"curve: "<<endl;
		//query_curve->print();

        if(flag == 1){
            break;
        }
        convert_2d_curve_to_vector(query_curve, grids[i], delta, hash_table_dimension,
		curve_dimension, &query_grid_curve, &query_item,max_coord);
        //cout <<"grid curve: "<<endl;
		//query_grid_curve->print();
		//cout <<"item: "<<endl;
		//query_item->print();
        P_value = hash_tables[i]->p(*(query_item->get_coordinates()),hash_table_dimension, table_size, w, k, bits_of_each_hash, M,  m_powers);
        int bucketes_checked=0;

        ret = hash_tables[i]->get_f_values_map()->equal_range(P_value);
		searched_items = 0;
    for (it = ret.first; it != ret.second; ++it) {
			if (searched_items >= M_f) {
				flag =1;
			}
      if(flag == 1){
        break;
            }
			if (check_for_identical_grid_flag == true) {
				if (it->second->get_corresponding_curve()->identical(query_curve) == false) {
					continue;
				}
			}

			double cur_distance = Curve_Grid_distance(query_curve, it->second->get_corresponding_curve());
			if (cur_distance < best_distance) {
				best = it->second->get_name();
				best_distance = cur_distance;
			}
			searched_items++;
		}
        unsigned nbuckets=hash_tables[i]->get_f_values_map()->bucket_count();
        for (unsigned y=0; y<nbuckets; y++) {
            if (searched_items >= M_f || bucketes_checked>=probes) {
				break;
			}
            for (auto it = hash_tables[i]->get_f_values_map()->begin(y);it!= hash_tables[i]->get_f_values_map()->end(y);it++){
	        	//std::cout << "[" << hash_table->get_f_values_map()->begin(i)->first << ":" << it->second->get_name() << "] ";
				bucket_value=hash_tables[i]->get_f_values_map()->begin(y)->first;
				break;
	    	}
            if (hammingDistance(P_value, bucket_value)==1){
				bucketes_checked++;
				ret = hash_tables[i]->get_f_values_map()->equal_range(bucket_value);
				for (it = ret.first; it != ret.second; ++it) {
					if (searched_items >= M_f) {
						break;
					}
					double cur_distance = Curve_Grid_distance((query_curve), (it->second->get_corresponding_curve()));//apostasi querry apo ta alla pou iparxoun sto bucket
					if (cur_distance < best_distance) {
						best = it->second->get_name();
						best_distance = cur_distance;
					}
					searched_items++;
				}
			}

        }

    }
    time = clock() - time;
    delete query_grid_curve;
    delete query_item;
	if (best != "") {
		query_result.set_best_distance(best_distance);
		query_result.set_time( ((double)time) / CLOCKS_PER_SEC );
		query_result.set_best_item(best);
	}
	else {
		query_result.set_best_distance(-1);
		query_result.set_time(-1);
		query_result.set_best_item("NULL");
	}

}





Hash_Table_Hypercube::Hash_Table_Hypercube(int table_size, int dimension, int w, int k){
	this->table_size=table_size;
    this->dimension = dimension;

	for (int i = 0; i < table_size; i++) {
		unordered_map<unsigned,int> *it=new unordered_map<unsigned,int>;
		g_value.push_back(it);
	}
    for (size_t y = 0; y < table_size; y++) {
        vector < vector <float>* > *it1 =new vector < vector <float>* >;
        for (int i = 0; i < k; i++){
            vector<float> *s = new vector<float>;
            random_float_vector(0, w, *s, dimension);
            it1->push_back(s);
        }
        s_array.push_back(it1);
    }
}

Hash_Table_Hypercube::~Hash_Table_Hypercube() {
	for (int i = 0; i < table_size; i++) {
		delete g_value[i];
	}

	for (vector < vector <float>* >* vv: s_array) {
		for (vector <float>* v : *vv) {
			delete v;
		}
		delete vv;
	}
}


void Hash_Table_Hypercube::insert(Curve *curve, unsigned g_value) {
	f_value.insert(pair<unsigned, Curve*>(g_value,curve) );
}

unsigned Hash_Table_Hypercube::p(vector<Type> x , int dimension, int table_size, int w, int k,
	int bits_of_each_hash, unsigned M, vector<unsigned>& m_powers) {
	int result=0;
	int p=0;
	//cout <<"total s_array size:"<<s_array.size()<<endl;
	//cout <<"table size"<<table_size<<endl;
	//getchar();
	for (int i = 0; i < table_size; i++) {
		//cout <<" p sa_ray size: "<<s_array[i]->size()<<endl;
		p = f_hash_function(x,dimension,w,k,bits_of_each_hash,M,m_powers,*s_array[i],g_value,i);
		result |= p << (table_size -i -1);
	}
	return result;
}
double Curve_Grid_hypercube::Curve_Grid_distance(Curve *curve1, Curve *curve2) {
	return DTW(curve1->get_points(), curve2->get_points());
}
