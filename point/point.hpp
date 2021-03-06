#ifndef POINT_H
#define POINT_H

#include <vector>

#include "../item/item.hpp"

using namespace std;

class Point {
public:
	Point() {}
	Point(Type x, Type y);

	void insert_coordinate(Type coordinate);
	vector<Type>& get_coordinates();
	void print_coordinates();
	Type get_coord(int i) const {return coordinates[i];}
	Type get_x() const {return coordinates[0];}
	Type get_y() const {return coordinates[1];}
	Type set_x(Type x) {coordinates.push_back(x);}
	Type set_y(Type y) {coordinates.push_back(y);}
	bool equals(Point *p);
private:
	vector<Type> coordinates;
};

#endif
