#ifndef TUPLE_H
#define TUPLE_H
#include <string>
using namespace std;
class Tuple {
public:
	Tuple() {}
	Tuple(int _x, int _y) :x(_x), y(_y){}
	int get_x() {return x;}
	int get_y() {return y;}
	int get_coord(int i) {return i == 0 ? x : y;}
	void set_x(int _x) {x = _x;}
	void set_y(int _y) {y = _y;}
private:
	int x, y;
};

#endif
