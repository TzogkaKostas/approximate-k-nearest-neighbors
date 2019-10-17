#ifndef CURVE_H
#define CURVE_H

#include "../point/point.hpp"

class Curve {
public:
	Curve(string name, vector<Point*> *points);
	Curve(vector<Point*> *points);
	~Curve();

	string get_name(){return name;}
	vector<Point*> get_points() {return *points;}
	void set_points(vector<Point*> *points) {this->points = points;}
	Curve* get_grid_curve() {return grid_curve;}
	void set_grid_curve(Curve *grid_curve) {this->grid_curve = grid_curve;}
	int get_length() {return points->size();}
	void print_curve();
	void print_grid_curve();
	void print_points();
	void insert_point(Point *point);
	bool identical(Curve *curve);

private:
	string name;
	vector<Point*> *points;
	Curve *grid_curve;
};

#endif