#include <vector>
#include <iostream>

using namespace std;

#include "../point/point.hpp"
#include "curve.hpp"

Curve::Curve(string name, vector<Point*>* points) {
	this->name = name;
	this->points = points;
	this->grid_curve = NULL;
}

Curve::Curve(vector<Point*>* points) {
	this->points = points;
	this->grid_curve = NULL;
}

Curve::~Curve() {
	for (size_t i = 0; i < points->size(); i++) {
		delete (*points)[i];
	}
	delete points;
	delete grid_curve;
}

void Curve::insert_point(Point *point) {
	points->push_back(point);
}

bool Curve::identical(Curve *curve) {
	typename vector<Point*>::iterator it1, it2;
	it1 = this->get_points().begin();
	it2 = curve->get_points().begin();

	for (; it1 != this->get_points().end() && it2 != curve->get_points().end() ; it1++, it2++) {
		if ((*it1)->equals(*it2) == false ) {
			return false;
		}
	}

	return true;
}

void Curve::print_curve() {
	cout <<name<<" ";
	cout <<"(";
	print_points();
	cout <<")"<<endl;
}

void Curve::print_grid_curve() {
	cout <<name<<" ";
	cout <<"(";
	grid_curve->print_points();
	cout <<")"<<endl;
}

void Curve::print_points() {
	for (Point *point : *points) {
		cout <<"(";
		point->print_coordinates();
		cout <<"), ";
	}
}