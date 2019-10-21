#ifndef HASH_H
#define HASH_H

#include <unordered_map>

#include "../curve/curve.hpp"

typedef Curve Content;

class Hash_Table {
public:
	Hash_Table(int dimension, int w, int k);
	~Hash_Table();

	//void insert_curve(Item *item, Curve *curve, int dimension, int w, int k, int bits_of_each_has,
	//		unsigned M, vector<unsigned> m_powers);
	void insert(Content* content, unsigned g_value);
	void print();
	unordered_multimap<unsigned, Content*>* get_map();
	vector<vector<float>*>& get_s_array();
private:
	unordered_multimap<unsigned, Content*> map;
	vector<vector<float>*> s_array;
};

#endif