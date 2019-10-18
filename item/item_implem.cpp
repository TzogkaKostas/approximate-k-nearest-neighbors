#include <vector>
#include <string>

using namespace std;

#include "../item/item.hpp"
#include "../helping_functions/helping_functions.hpp"

Item::Item(string name, vector<Type>* coordinates) {
	this->name = name;
	this->coordinates = coordinates;
}

Item::Item(vector<Type>* coordinates) {
	this->coordinates = coordinates;
}

Item::~Item() {
	delete coordinates;
}

void Item::print() { 
	print_vector(*coordinates);
}