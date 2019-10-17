#include <vector>
#include <iostream>

using namespace std;

#include "point.hpp"

Point::Point(vector<Type> *coordinates) {
	this->coordinates = coordinates;
}
Point::Point(Type x, Type y) {
	this->coordinates = new vector<Type>;
	this->coordinates->push_back(x);
	this->coordinates->push_back(y);
}

Point::Point(const Point &p) {
	this->coordinates = p.get_coordinates();
}

Point::~Point() {
	delete coordinates;
}

vector<Type>* Point::get_coordinates() const {
	return coordinates;
}

void Point::set_coordinates(vector<Type> *coordinates) {
	this->coordinates = coordinates;
}

bool Point::equals(Point *p) {
	return get_x() == p->get_x() && get_y() == p->get_y() ? true : false;
}

void Point::print_coordinates() {
	for (size_t i = 0; i < coordinates->size(); i++) {
		cout <<(*coordinates)[i]<<", ";
	}
}