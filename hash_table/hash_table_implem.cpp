#include <vector>
#include <string>
#include <unordered_map>

using namespace std;

#include "../item/item.hpp"
#include "../helping_functions/helping_functions.hpp"
#include "../hash_table/hash_table.hpp"

Hash_Table::Hash_Table(int dimension, int w, int k) {
	for (size_t i = 0; i < k; i++){
		vector<float> *s = new vector<float>;
		random_float_vector(0, w, *s, dimension);
		s_array.push_back(s);
	}
}

Hash_Table::~Hash_Table() {
	for (size_t i = 0; i < s_array.size(); i++){
		delete s_array[i];
	}
}

//void Hash_Table::insert_curve(Item *item, Curve *curve, int dimension, int w, int k,
//			int bits_of_each_hash, unsigned M, vector<unsigned> m_powers) {
//	unsigned g_value, position;
//	g_value = g_hash_function(*(item->get_coordinates()), dimension,
//			w, k, bits_of_each_hash, M, s_array, m_powers);
//
//	map.insert({g_value, curve});
//}

void Hash_Table::insert(Content *content, unsigned g_value) {
	map.insert({g_value, content});
}

unordered_multimap<unsigned, Content*>* Hash_Table::get_map() {
	return &map;
}
vector<vector<float>*>& Hash_Table::get_s_array() {
	return s_array;
}

void Hash_Table::print() {
	//for (size_t i = 0; i < get_table_size(); i++) {
	//	if (buckets[i].empty() == false) {
	//		cout <<"Bucket:" << i <<endl;
	//		print_list(buckets[i]);
	//	}
	//}   
}