#ifndef TUPLE_H
#define TUPLE_H

class Tuple {
public:
	Tuple() {}
	Tuple(int _x, int _y) :x(_x), y(_y){}
	int get_x() {return x;}
	int get_y() {return y;}
	void set_x(int _x) {x = _x;}
	void set_y(int _y) {y = _y;}
private:
	int x, y;
};

#endif
