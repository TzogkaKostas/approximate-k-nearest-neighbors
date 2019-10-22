#ifndef ITEM_H
#define ITEM_H

#include <vector>
#include <string>
using namespace std;

typedef float Type;

class Item {
public:
	Item() {}
	Item(vector<Type>* coordinates);
	Item(string name, vector<Type>* coordinates);
	~Item();

	string get_name() const {return name;}
	void set_name(string name) {this->name = name;}
	vector<Type>* get_coordinates() {return coordinates;}
	void set_coordinates(vector<Type>* coordinates) {this->coordinates = coordinates;}

	void print();
private:
	string name;
	vector<Type> *coordinates;
};

#endif
