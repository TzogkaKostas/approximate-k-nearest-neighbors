#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
#include <list>
#include <algorithm>
#include <cstring>
#include <vector>
#include <string>
#include <list>
#include <vector>

using namespace std;

#include "../item/item.hpp"
#include "../lsh/lsh.hpp"
#include "../query_result/query_result.hpp"
#include "../hash_table/hash_table.hpp"
#include "../query_result/query_result.hpp"
#include "../helping_functions/helping_functions.hpp"

void read_command_line_arguments(char *argv[], int& argc, string& input_file, string& query_file,
	string& output_file, int& k, int& L, int& w, int& st, bool& check_for_identical_grid_flag) {
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-d") == 0) {
			input_file = argv[i + 1];
		}
		else if (strcmp(argv[i], "-q") == 0) {
			query_file = argv[i + 1];
		}
		else if (strcmp(argv[i], "-o") == 0) {
			output_file = argv[i + 1];
		}
		else if (strcmp(argv[i], "-k") == 0) {
			k = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-L") == 0) {
			L = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-w") == 0) {
			w = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-st") == 0) {
			st = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "--identical") == 0) {
			check_for_identical_grid_flag = atoi(argv[i + 1]);
		}
	}
}

void read_command_line_arguments(char *argv[], int& argc, string& input_file, string& query_file,
	string& output_file, int& k, int& L, int& w, int& st) {
	bool x;
	//read_command_line_arguments(argv, argc, input_file, query_file, output_file,
	//		k, L, w, st, &x);
}

void delete_curves(list<Curve*> curves) {
	for (Curve *curve : curves) {
		delete curve;
	}
}

int read_2d_curves_from_file(string file_name, list<Curve*>& curves, int& max_length) {
	string line, coordinate, name, tuple, part1, part2;
	int length;
	Type x, y;
	vector<Type> *coordinates;
	Point *point;
	Curve *curve;
	vector<Point*> *points;

	ifstream inputfile(file_name.c_str());
	if (inputfile.is_open() == false) {
		cout <<"File opening error: "<<file_name<<endl;
		return 1;
	}
	else {
		max_length = -1;
		while(getline(inputfile, line)) {
			istringstream iss (line);
			iss >> name;
			iss >> length;
			max_length = max(max_length, length);

			points = new vector<Point*>;
			while ( iss >> part1 >> part2) {
				tuple = part1 + part2;
				sscanf(tuple.c_str(), "(%f, %f)", &x, &y);
				coordinates = new vector<Type>;
				coordinates->push_back(x);
				coordinates->push_back(y);
				point = new Point(coordinates);
				points->push_back(point);
			}
			curve = new Curve(name, points);
			curves.push_back(curve);
		}
		inputfile.close();
	}
	return 0;
}

void print_curves(list<Curve*> curves)  {
	for (Curve *curve : curves) {
		curve->print_curve();
		cout <<endl;
	}
}

Type DTW(vector<Point*> p, vector<Point*> q) {
	int m1 = p.size();
	int m2 = q.size();

	Type dtw_array[m1][m2];

	//initialize (0,0) point of the array
	dtw_array[0][0] = manhattan_distance_2d(p[0], q[0]);

	//initialize first row of the array
	for (size_t j = 1; j < m2; j++) {
		dtw_array[0][j] = manhattan_distance_2d( p[0], q[j] ) + dtw_array[0][j - 1];
	}

	//initialize first column of the array
	for (size_t i = 1; i < m1; i++) {
		dtw_array[i][0] = manhattan_distance_2d( p[i], q[0] ) + dtw_array[i - 1][0];
	}

	//calculate the rest cells of the array( (1,1) is calculated twice, but that's ok )
	for (size_t i = 1; i < m1; i++) {
		for (size_t j = 1; j < m2; j++) {
			dtw_array[i][j] = manhattan_distance_2d( p[i], q[j] ) +
				min({dtw_array[i - 1][j - 1], dtw_array[i][j - 1], dtw_array[i - 1][j]});
		}
	}
	return dtw_array[m1 - 1][m2 - 2];
}

float manhattan_distance_2d(Point *p, Point *q) {
	return abs(p->get_x() - q->get_x()) + abs(p->get_y() - q->get_y());
}

unsigned g_hash_function(vector<Type> x , int dimension, int w, int k, 
	int bits_of_each_hash, unsigned M, vector<vector<float>*>& s_array, vector<unsigned>& m_powers) {		

	unsigned hash_value, total_hash_value = 0;
	for (int i = 0; i < k; ++i) {
		hash_value = 0;
		hash_value = hash_function(x, dimension, w, M, *(s_array[i]), m_powers);
		total_hash_value |= hash_value << (32 - (i + 1)*bits_of_each_hash);
	}
	return total_hash_value;
}

unsigned hash_function(vector<Type> x, int dimension, int w, unsigned M,
	vector<float>& s, vector<unsigned>& m_powers) {
	unsigned sum = 0;
	vector<int> a;
	for (size_t i = 0; i < dimension; i++) {
		a.push_back( floor( (x[i] - s[i]) / w) );
	}

	sum = mod( mul_mod(a[0], m_powers[dimension - 1 - 0], M), M);
	for (size_t i = 1; i < a.size(); i++) {
		sum = add_mod( mul_mod(a[i], m_powers[dimension - 1 - i], M), sum, M);
	}

	return sum;
}


void print_ann_results(Query_Result ann_result) {
	if (ann_result.get_name() == "NULL") {
		cout <<"NOT FOUND"<<endl;
		return;
	}
	cout <<"Approx Nearest neighbor:"<<ann_result.get_name()<<endl;
	cout <<"distanceLSH:"<<ann_result.get_best_distance()<<endl;
	cout <<"tLSH:"<<ann_result.get_time()<<endl;
}

void print_exhaustive_search_results(Query_Result exhaustive_result) {
	if (exhaustive_result.get_name() == "NULL") {
		cout <<"NOT FOUND"<<endl;
		return;
	}
	cout <<"Nearest neighbor:"<<exhaustive_result.get_name()<<endl;
	cout <<"distanceTrue:"<<exhaustive_result.get_best_distance()<<endl;
	cout <<"tTrue:"<<exhaustive_result.get_time()<<endl;
}

void print_parameters(int L, int k, int w, int search_threshold, int dimension) {
	cout <<"L: "<<L<<endl;
	cout <<"k: "<<k<<endl;
	cout <<"w: "<<w<<endl;
	cout <<"search_threshold: "<<search_threshold<<endl;
	cout <<"dimension: "<<dimension<<endl<<endl;;
}

void get_vector_from_line(string line, Item& item) {
	string name;
	Type coordinate;

	istringstream iss (line);
	iss >> name;
	vector<Type> *coordinates = new vector<Type>;
	while ( iss >> coordinate ) {
		coordinates->push_back(coordinate);
	}
	item.set_name(name);
	item.set_coordinates(coordinates);
}

void read_vectors_from_file(string file_name, list<Item*>& items) {
	string line, name;
	Type coordinate;

	ifstream inputfile(file_name.c_str());
	if (inputfile.is_open() == false) {
		cout <<"File opening error: "<<file_name<<endl;
		return;
	}
	else {
		while(getline(inputfile, line)) {
			istringstream iss (line);
			iss >> name;
			vector<Type> *coordinates = new vector<Type>;
			while ( iss >> coordinate ) {
				coordinates->push_back(coordinate);
			}
			Item *item = new Item(name, coordinates);
			items.push_back(item);
		}
		inputfile.close();
	}
	return;
}

unsigned long long int manhattan_distance(vector<Type> x1, vector<Type> x2) {

	unsigned long long int sum = 0;
	for (size_t i = 0; i < x1.size(); i++) {
		sum += abs(x1[i] - x2[i]);
	}
	
	return sum;
}


void exhaustive_curve_search(list<Curve*> *curves, Curve *query, Query_Result& query_result) {
	unsigned best_distance = numeric_limits<unsigned>::max();
	string best = "";

	time_t time;
	time = clock();
	for (Curve *curve : *curves) {
		unsigned cur_distance = DTW(curve->get_points(), query->get_points());
		if (cur_distance < best_distance) {
			best = curve->get_name();
			best_distance = cur_distance;
		}
	}
	time = clock() - time;
	
	if (best != "") {
		query_result.set_best_distance(best_distance);
		query_result.set_time( ((double)time) / CLOCKS_PER_SEC );
		query_result.set_best_item(best);
	}
	else {
		query_result.set_best_distance(-1);
		query_result.set_time(-1);
		query_result.set_best_item("NULL");
	}
}


void exhaustive_search(list<Item*> *items, Item *query, Query_Result& query_result) {

	unsigned best_distance = numeric_limits<unsigned>::max();
	string best = "-";

	time_t time;
	time = clock();
	for (Item *item : *items) {
		unsigned cur_distance = manhattan_distance(*(query->get_coordinates()), *(item->get_coordinates()) );
		if (cur_distance < best_distance) {
			best = item->get_name();
			best_distance = cur_distance;
		}
	}
	time = clock() - time;
	
	if (best != "-") {
		query_result.set_best_distance(best_distance);
		query_result.set_time( ((double)time) / CLOCKS_PER_SEC );
		query_result.set_best_item(best);
	}
	else {
		query_result.set_best_distance(-1);
		query_result.set_time(-1);
		query_result.set_best_item("NULL");
	}
}

void delete_items(list<Item*> items) {
	for (Item *item : items) {
		delete item;
	}
}


unsigned mod(long long a, long long b) {
	return (a%b + b)%b;
}

unsigned add_mod(long long a, long long b, long long m) {
	if ( 0 == b ) return a;

	b = m - b;
	if ( a>=b ) {
		return a - b;
	}
	else {
		return m - b + a;
	}
}

uint32_t mul_mod(long long a, long long b, long long m) { 
	long double x;
	uint32_t c;
	int32_t r;
	if (a >= m) a %= m;
	if (b >= m) b %= m;

	x = a;
	c = x * b / m;
	r = (int64_t)(a * b - c * m) % (int64_t)m;
	return r < 0 ? r + m : r;
}

uint32_t mul_mod2(uint32_t a, uint32_t b, uint32_t m) { 
	long double x;
	uint32_t c;
	int32_t r;
	if (a >= m) a %= m;
	if (b >= m) b %= m;
	x = a;
	c = x * b / m;
	r = (int64_t)(a * b - c * m) % (int64_t)m;
	return r < 0 ? r + m : r;
}

uint32_t pow_mod(uint32_t a, uint32_t b, uint32_t m) { 
	uint32_t r = m==1?0:1;
	while (b > 0) {
		if (b & 1)
			r = mul_mod2(r, a, m);
		b = b >> 1;
		a = mul_mod2(a, a, m);
	}
	return r;
}

void random_float_list(double from, double to, list<float>& random_list, int size) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<float> dis(from, to);
	for (size_t i = 0; i < size; i++) {
		random_list.push_back(dis(gen));
	}
}

void random_float_vector(double from, double to, vector<float>& random_vector, int size) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<float> dis(from, to);
	for (size_t i = 0; i < size; i++) {
		random_vector.push_back(dis(gen));
	}
}


int max_from_list(list<int> my_list) {
	return *max_element(my_list.begin(), my_list.end());
}

int max_from_list(list<float> my_list) {
	return *max_element(my_list.begin(), my_list.end());
}

int max_from_vector(vector<float> my_vector) {
	return *max_element(my_vector.begin(), my_vector.end());
}

int max_from_vector(vector<int> my_vector) {
	return *max_element(my_vector.begin(), my_vector.end());
}


int find_dimension_from_file(string file_name) {
	string line, word;
	int counter;

	ifstream inputfile(file_name.c_str());
	if (inputfile.is_open() == false) {
		cout <<"File opening error: "<<file_name<<endl;
		return 1;
	}
	else {
		getline(inputfile, line);
		istringstream iss (line);
		counter = 0;
		while ( iss >> word ) {
			counter++;
		}
		inputfile.close();
	}
	return counter - 1;
}

int get_lines_of_file(string file_name) {
	ifstream inFile(file_name); 
 	return count(istreambuf_iterator<char>(inFile), istreambuf_iterator<char>(), '\n');
}


int print_input_file(string file_name) {
	string line, coordinate, name;

	ifstream inputfile(file_name.c_str());
	if (inputfile.is_open() == false) {
		cout <<"File opening error: "<<file_name<<endl;
		return 1;
	}
	else {
		while(getline(inputfile, line)) {
			istringstream iss (line);
			iss >> name;
			cout << name << " ";
			while ( iss >> coordinate ) {
				cout << coordinate << " ";
			}
			cout <<endl;
		}
		inputfile.close();
	}
	return 0;
}

int max_from_array(int array[], int size) {

	int max = array[0];
	for (size_t i = 1; i < size; i++) {
		  if (array[i] > max) {
			max  = array[i];
		}
	}
	return max;
}

int max_from_array(unsigned array[], int size) {

	unsigned max = array[0];
	for (size_t i = 1; i < size; i++) {
		  if (array[i] > max) {
			max  = array[i];
		}
	}
	return max;
}

float max_from_array(float array[], int size) {

	float max = array[0];
	for (size_t i = 1; i < size; i++) {
		  if (array[i] > max) {
			max  = array[i];
		}
	}
	return max;
}

double max_from_array(double array[], int size) {

	double max = array[0];
	for (size_t i = 1; i < size; i++) {
		  if (array[i] > max) {
			max  = array[i];
		}
	}
	return max;
}

void random_double_array(double from, double to, double *random_array, int size) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(from, to);
	for (size_t i = 0; i < size; i++) {
		random_array[i] =  dis(gen);
	}
}

void random_float_array(double from, double to, float *random_array, int size) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<float> dis(from, to);
	for (size_t i = 0; i < size; i++) {
		random_array[i] =  dis(gen);
	}
}

void random_int_list(double from, double to, list<int>& random_list, int size) {
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dis(from, to);
	for (size_t i = 0; i < size; i++) {
		random_list.push_back(dis(gen));
	}
}

double random_double(double from, double to) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(from, to);
	return dis(gen);
}

float random_float(double from, double to) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<float> dis(from, to);
	return dis(gen);
}

int random_int(int from, int to) {
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dis(from, to);
	return dis(gen);
}

void print_list(list<int> my_list) {
	for (int i : my_list) {
		cout << i << " ";
	}
	cout <<endl;
}

void print_vector(vector<int> my_list) {
	for (int i : my_list) {
		cout << i << " ";
	}
	cout <<endl;
}

void print_vector(vector<float> my_list) {
	for (float i : my_list) {
		cout << i << " ";
	}
	cout <<endl;
}


void print_list(list<float> my_list) {
	for (float i : my_list) {
		cout << i << " ";
	}
	cout <<endl;
}

void print_array(double array[], int size) {
	for (size_t i = 0; i < size; i++) {
		cout << array[i] << " ";
	}
	cout <<endl;
}

void print_array(float array[], int size) {
	for (size_t i = 0; i < size; i++) {
		cout << array[i] << " ";
	}
	cout <<endl;
}

void print_array(int array[], int size) {
	for (size_t i = 0; i < size; i++) {
		cout << array[i] << " ";
	}
	cout <<endl;
}