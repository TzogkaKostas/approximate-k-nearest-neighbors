#ifndef RES_H
#define RES_H

#include <string>
using namespace std;

class Query_Result {
public:
	Query_Result() :best_item("-"){}
	string get_name() {return best_item;}
	long long int get_best_distance() {return best_distance;}
	double get_time() {return time;}

	void set_best_item(string name) {this->best_item = name;}
	void set_best_distance(long long int best_distance) {this->best_distance = best_distance;}
	void set_time(double time) {this->time = time;}
private:
	string best_item;
	long long int best_distance;
	double time;
};

#endif