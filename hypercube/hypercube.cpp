#include "hypercube.hpp"


using namespace std;

typedef float Type;


Hash_Table::Hash_Table(int table_size, int dimension, int w, int k){
	this->table_size=table_size;
	for (int i = 0; i < k; i++){
		vector<float> *s = new vector<float>;
		random_float_vector(0, w, *s, dimension);
		s_array.push_back(s);
	}
	for (int i = 0; i < table_size; i++) {
		unordered_map<unsigned,int> *it=new unordered_map<unsigned,int>;
		g_value.push_back(it);
	}

}

Hash_Table::~Hash_Table() {
	for (int i = 0; i < table_size; i++) {
		delete g_value[i];
	}

}


void Hash_Table::insert(Item*item,int dimension, int w, int k,int bits_of_each_hash, unsigned M,vector<unsigned> m_powers) {
	unsigned index =p(*item->get_coordinates(),dimension,table_size,w,k,bits_of_each_hash,M,m_powers);
	f_value.insert(pair<unsigned, Item*>(index,item) );

}

unsigned Hash_Table::p(vector<Type> x , int dimension, int table_size, int w, int k,
	int bits_of_each_hash, unsigned M, vector<unsigned>& m_powers) {
	int result=0;
	int p=0;
	for (int i = 0; i < table_size; i++) {

		p = f_hash_function(x,dimension,w,k,bits_of_each_hash,M,m_powers,s_array,g_value,i);
		result = result ^ p;
		result = result <<1;
	}
	return result;
}



void Hypercube::insert_item(Item *item){
		hash_table->insert(item,this->dimension,this->w,this->k,-this->bits_of_each_hash,this->M,this->m_powers);
}

Hypercube::Hypercube(int hash_table_size, int dimension, int w, int k,unsigned m,unsigned M){
	this->dimension=dimension;
	this->w=w;
	this->k=k;
	this->bits_of_each_hash=32/k;
	this->m=m;
	this->table_size=hash_table_size;
	this->M_f=M;
	this->hash_table = new Hash_Table(hash_table_size, dimension, w, k);
	if (k == 1) {
		this->M = numeric_limits<unsigned>::max();
	}
	else {
		this->M = pow(2, bits_of_each_hash);
	}
	for (int i = 0; i < dimension; i++) {
		m_powers.push_back( pow_mod(m, i, M) );
	}

}

void Hypercube::ANN(Item *query, unsigned prompt, Query_Result& query_result){
	int searched_items;
	unsigned best_distance = numeric_limits<unsigned>::max();
	unsigned  F_value;
	string best = "";
	unordered_multimap<unsigned, Item*> *map;
	pair <unordered_multimap<unsigned, Item*>::iterator, unordered_multimap<unsigned,Item*>::iterator> ret;
	unordered_multimap<unsigned, Item*>::iterator it;

	int bucketes_checked=0;
	time_t time;
	time = clock();
	F_value = hash_table->p(*(query->get_coordinates()),dimension,
	                        table_size, w, k, bits_of_each_hash, M,  m_powers);

	map = hash_table->get_f_values_map();
	ret = map->equal_range(F_value);
	searched_items = 0;

	unsigned bucket_value;

		for (it = ret.first; it != ret.second; ++it) {
			if (searched_items >= M_f) {
				break;
			}
			unsigned cur_distance = Hypercube_distance(query, it->second);//apostasi querry apo ta alla pou iparxoun sto bucket
			if (cur_distance < best_distance) {
				best = it->second->get_name();
				best_distance = cur_distance;
			}
			searched_items++;
		}
		unsigned nbuckets=hash_table->get_f_values_map()->bucket_count();
		for (unsigned i=0; i<nbuckets; i++) {
			if (searched_items >= M_f || searched_items>=prompt) {
				break;
			}

	  	  	for (auto it = hash_table->get_f_values_map()->begin(i);it!= hash_table->get_f_values_map()->end(i);it++){
	        	//std::cout << "[" << hash_table->get_f_values_map()->begin(i)->first << ":" << it->second->get_name() << "] ";
				bucket_value=hash_table->get_f_values_map()->begin(i)->first;
				break;
	    	}
			if (hammingDistance(F_value, bucket_value)==1){
				bucketes_checked++;
				ret = map->equal_range(bucket_value);
				for (it = ret.first; it != ret.second; ++it) {
					if (searched_items >= M_f) {
						break;
					}
					unsigned cur_distance = Hypercube_distance(query, it->second);//apostasi querry apo ta alla pou iparxoun sto bucket
					if (cur_distance < best_distance) {
						best = it->second->get_name();
						best_distance = cur_distance;
					}
					searched_items++;
				}
			}
	    }


	time = clock() - time;

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
void Hypercube::range_search(Item *query, unsigned prompt, float radious,
		list<Item*>& range_items, Query_Result& query_result) {
		int searched_items;
		unsigned best_distance = numeric_limits<unsigned>::max();
		unsigned  F_value;
		string best = "";
		unordered_multimap<unsigned, Item*> *map;
		pair <unordered_multimap<unsigned, Item*>::iterator, unordered_multimap<unsigned,Item*>::iterator> ret;
		unordered_multimap<unsigned, Item*>::iterator it;

		int bucketes_checked=0;
		time_t time;
		time = clock();
		F_value = hash_table->p(*(query->get_coordinates()),dimension,
		                        table_size, w, k, bits_of_each_hash, M,  m_powers);

		map = hash_table->get_f_values_map();
		ret = map->equal_range(F_value);
		searched_items = 0;

		unsigned bucket_value;

			for (it = ret.first; it != ret.second; ++it) {
				if (searched_items >= M_f) {
					break;
				}
				unsigned cur_distance = Hypercube_distance(query, it->second);//apostasi querry apo ta alla pou iparxoun sto bucket
				if (cur_distance < radious) {
					range_items.push_back(it->second);
				}
				searched_items++;
			}
			unsigned nbuckets=hash_table->get_f_values_map()->bucket_count();
			for (unsigned i=0; i<nbuckets; i++) {
				if (searched_items >= M_f || searched_items>=prompt) {
					break;
				}

		  	  	for (auto it = hash_table->get_f_values_map()->begin(i);it!= hash_table->get_f_values_map()->end(i);it++){
		        	//std::cout << "[" << hash_table->get_f_values_map()->begin(i)->first << ":" << it->second->get_name() << "] ";
					bucket_value=hash_table->get_f_values_map()->begin(i)->first;
					break;
		    	}
				if (hammingDistance(F_value, bucket_value)==1){
					bucketes_checked++;
					ret = map->equal_range(bucket_value);
					for (it = ret.first; it != ret.second; ++it) {
						if (searched_items >= M_f) {
							break;
						}
						unsigned cur_distance = Hypercube_distance(query, it->second);//apostasi querry apo ta alla pou iparxoun sto bucket
						if (cur_distance < radious) {
							range_items.push_back(it->second);
						}
						searched_items++;
					}
				}
		    }


		time = clock() - time;

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

unsigned long long int Hypercube::Hypercube_distance(Item *item1, Item *item2) {
	return manhattan_distance(*(item1->get_coordinates()), *(item2->get_coordinates()));
}
