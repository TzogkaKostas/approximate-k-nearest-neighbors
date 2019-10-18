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
	Curve* get_corresponding_curve() {return corresponding_curve;}
	void set_corresponding_curve(Curve *corresponding_curve) {this->corresponding_curve = corresponding_curve;}
	int get_length() {return points->size();}
	void print();
	void print_corresponding_curve();
	void print_points();
	void insert_point(Point *point);
	bool identical(Curve *curve);

private:
	string name;
	vector<Point*> *points;
	Curve *corresponding_curve;
};

#endif