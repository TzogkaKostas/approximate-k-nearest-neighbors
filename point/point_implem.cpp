#include <vector>
#include <iostream>

using namespace std;

#include "point.hpp"


Point::Point(Type x, Type y) {
	this->coordinates.push_back(x);
	this->coordinates.push_back(y);
}
vector<Type>& Point::get_coordinates()  {
	return coordinates;
}

void Point::insert_coordinate(Type coordinate) {
	this->coordinates.push_back(coordinate);

}

bool Point::equals(Point *p) {
	return get_x() == p->get_x() && get_y() == p->get_y() ? true : false;
}

void Point::print_coordinates() {
	for (size_t i = 0; i < coordinates.size(); i++) {
		cout <<coordinates[i]<<", ";
	}
}