#ifndef POINT_H
#define POINT_H

#include <vector>

#include "../item/item.hpp"

class Point {
public:
	Point(vector<Type> *coordinates);
	Point(Type x, Type y);
	Point(const Point &p);
	~Point();

	vector<Type>* get_coordinates() const;
	void set_coordinates(vector<Type> *coordinates);
	void print_coordinates();
	Type get_x() const {return (*coordinates)[0];}
	Type get_y() const {return (*coordinates)[1];}
	Type set_x(Type x) {coordinates->push_back(x);}
	Type set_y(Type y) {coordinates->push_back(y);}
	bool equals(Point *p);
private:
	vector<Type> *coordinates;
};

#endif