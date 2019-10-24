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
#include "../Tuple/tuple.hpp"

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_int_distribution<> dis(0, 1);

void convert_2d_curve_to_vector_by_projection(vector<Tuple*>& U, float **G_matrix, Curve *curve,
		int G_rows, int G_cols, Item **item) {
	matrix_multiplication(U, G_matrix, curve, G_rows, G_cols, item);
}

void matrix_multiplication(vector<Tuple*>& U, float **G_matrix, Curve *curve,
		int G_rows, int G_cols, Item **item) {

	float sum;
	int position_of_curve;
	vector<Type> *results_points;

	//cout <<"Grows: "<<G_rows<<endl;
	//cout <<"G_cols: "<<G_cols<<endl;
	//cout <<"u.size: "<<U.size()<<endl;
	//curve->print();
	//for(Tuple *t : U) {
	//	cout <<"("<<t->get_x()<<", "<<t->get_y()<<") ";
	//}
	//cout <<endl;

	//for (size_t i = 0; i < G_rows; i++) {
	//	for (size_t j = 0; j < G_cols; j++) {
	//		G_matrix[i][j] = 0.5;
	//		cout <<G_matrix[i][j]<<" ";
	//	}
	//	cout << endl;
	//}

	vector<Point*>points = curve->get_points();
	results_points = new vector<Type>;
	for (size_t U_i = 0; U_i < U.size(); U_i++) {
		for(int i = 0; i < G_rows; ++i) { //for every row of G
        	for(int j = 0; j < 1; ++j) { //for every column of U
				sum = 0.0;
            	for(int k = 0; k < G_cols; ++k) {
					position_of_curve = U[U_i]->get_x();
					sum += G_matrix[i][k] * points[position_of_curve]->get_coord(k);
				}
				results_points->push_back(sum);
            }
		}
	}

	//for (Type p : *results_points) {
	//	cout <<p<<" ";
	//}
	//cout <<endl;
	//cout <<"done"<<endl;

	*item = new Item(results_points);
}

void convert_2d_curve_to_vector(Curve *curve, Point *t, int delta, int dimension,
		int curve_dimension, Curve **grid_curve, Item **item) {
	snap_curve(curve, t, grid_curve, delta);
	fill_curve(*grid_curve, dimension/curve_dimension - (*grid_curve)->get_length()); //padding
	(*grid_curve)->set_corresponding_curve(curve);
	zip_points(*grid_curve, item);
}

void read_command_line_arguments_hypercube(char *argv[], int& argc,string& input_file, string& query_file,
	string& output_file, int& k, int& M,int& probes,int &flag){
	int opt;
	while((opt = getopt(argc, argv, "d:q:o:k:M:p:")) != -1)
    {
      switch(opt){
          case 'd':
            	input_file= optarg;
          break;
          case 'q':
            	query_file=optarg;
          break;
		  case 'o':
		  			output_file=optarg;
		  break;
          case 'k':
               	k=atoi(optarg);
				flag=1;
          break;
          case 'M':{
		  		M=atoi(optarg);
		  }
          break;
          case 'p':
		  		probes=atoi(optarg);
          break;
      }
  }

}

void read_command_line_arguments_hypercube_grid(char *argv[], int& argc,string& input_file, string& query_file,
	string& output_file, int& k, int& M,int& probes,int &L,int &flag ){
	int opt;
	while((opt = getopt(argc, argv, "d:q:o:k:M:p:L:")) != -1)
    {
      switch(opt){
          case 'd':
            	input_file= optarg;
          break;
          case 'q':
            	query_file=optarg;
          break;
		  case 'o':
		  		output_file=optarg;
		  break;
          case 'k':
               	k=atoi(optarg);
				flag=1;
          break;
          case 'M':{
		  		M=atoi(optarg);
		  }
          break;
          case 'p':
		  		probes=atoi(optarg);
          break;
		  case 'L':
		  		L=atoi(optarg);
          break;
      }
  }

}
void read_command_line_arguments_hypercube_projection(char* argv[],int &argc,string &input_file,string &query_file,int& k, int& M,int & probes,float &e,string &output_file,int &flag ){
	int opt;
	while((opt = getopt(argc, argv, "d:q:o:k:M:p:e:")) != -1)
    {
      switch(opt){
          case 'd':
            	input_file= optarg;
          break;
          case 'q':
            	query_file=optarg;
          break;
		  case 'o':
		  		output_file=optarg;
		  break;
          case 'k':
               	k=atoi(optarg);
				flag=1;
          break;
          case 'M':{
		  		M=atoi(optarg);
		  }
          break;
          case 'p':
		  		probes=atoi(optarg);
          break;
		  case 'e':
		  		e=atof(optarg);
          break;
      }
  }

}

unsigned f_hash_function(vector<Type> x , int dimension,int w, int k,
	int bits_of_each_hash, unsigned M, vector<unsigned>& m_powers,vector<vector<float>*>& s_array,vector< unordered_map<unsigned,int> *>&g_value,int f){
	int y;
	srand(time(0));
	unsigned g;
	g=0;
	g=g_hash_function(x , dimension, w,  k,bits_of_each_hash,  M,s_array,m_powers);
	//cout << g_value[f]->find(g)->first<<endl;
	if(g_value[f]->find(g) == g_value[f]->end()){// not found
		y=dis(gen);
		g_value[f]->insert(make_pair(g,y));
		return y;
	}
	else{
		return g_value[f]->find(g)->second;
	}


}

void snap_curve(Curve *curve, Point *t, Curve **grid_curve, int delta) {

	Point *snapped_point = NULL;

	vector<Point*> *snapped_points = new vector<Point*>;
	for(Point *point : curve->get_points() ) {
		get_snapped_point(point, delta, t, &snapped_point);
		//cout <<"point:"<<endl;
		//point->print_coordinates();
		//cout <<endl;
		//cout <<"snapped point:"<<endl;
		//snapped_point->print_coordinates();
		//cout <<endl;
		//getchar();

		//consecutive duplicate
		if (snapped_points->size() > 0) {
			if (snapped_points->back()->equals(snapped_point) == true) {
				delete snapped_point;
				continue;
			}
		}
		snapped_points->push_back(snapped_point);
	}
	*grid_curve = new Curve(curve->get_name(), snapped_points);

}

void get_snapped_point(Point *point, int delta, Point *t, Point **snapped_point) {
	/*
			one cell, shifted by t, of the big Grid
		up_left -------------- up_right
				|\ d1	 d2 /|
				| \		   / |
				|	 point   |
				|  / d3	  \d4|
		down_left ------------- down_left
	*/
	//cout <<"t_x: "<<t->get_x()<<endl;
	//cout <<"t_y: "<<t->get_y()<<endl;
	//cout <<"point_x: "<<point->get_x()/delta<<endl;
	//cout <<"point_y: "<<point->get_y()/delta<<endl;
	//cout <<"delta: "<<delta<<endl;


	float up_left_x = floor( point->get_x()/delta )*delta + t->get_x();
	float up_left_y = ceil( point->get_y()/delta )*delta + t->get_y();
	//cout <<"up_left_x: "<<up_left_x<<endl;
	//cout <<"up_left_y: "<<up_left_y<<endl;

	float up_right_x = ceil( point->get_x()/delta)*delta + t->get_x();
	float up_right_y = ceil( point->get_y()/delta)*delta + t->get_y();
	//cout <<"up_right_x: "<<up_right_x<<endl;
	//cout <<"up_right_y: "<<up_right_y<<endl;

	float down_left_x = floor( point->get_x()/delta)*delta + t->get_x();
	float down_left_y = floor( point->get_y()/delta)*delta + t->get_y();

	//cout <<"down_left_x: "<<down_left_x<<endl;
	//cout <<"down_left_y: "<<down_left_y<<endl;

	float down_right_x = ceil( point->get_x()/delta)*delta + t->get_x();
	float down_right_y = floor( point->get_y()/delta)*delta + t->get_y();

	//cout <<"down_right_x: "<<down_right_x<<endl;
	//cout <<"down_right_y: "<<down_right_y<<endl;

	//maxhattan distances
	float d1 = abs(up_left_x - point->get_x()) + abs(up_left_y - point->get_y());
	float d2 = abs(up_right_x - point->get_x()) + abs(up_right_y - point->get_y());
	float d3 = abs(down_left_x - point->get_x()) + abs(down_left_y - point->get_y());
	float d4 = abs(down_right_x - point->get_x()) + abs(down_right_y - point->get_y());

	float min_d = d1;
	float best_x = up_left_x;
	float best_y = up_left_y;

	if (d2 < min_d) {
		min_d = d2;
		best_x = up_right_x;
		best_y = up_right_y;
	}

	if (d3 < min_d) {
		min_d = d3;
		best_x = down_left_x;
		best_y = down_left_y;
	}

	if (d4 < min_d) {
		min_d = d4;
		best_x = down_right_x;
		best_y = down_right_y;
	}
	*snapped_point = new Point(best_x, best_y);
}

void fill_curve(Curve *curve, int pad_length) {
	for (size_t i = 0; i < pad_length; i++) {
		Point *point = new Point(65536, 65536);
		curve->insert_point(point);
	}
}

void zip_points(Curve *grid_curve, Item **item) {
	vector<Type> *coordinates = new vector<Type>;
	for(Point *point : grid_curve->get_points() ) {
		coordinates->push_back(point->get_x());
		coordinates->push_back(point->get_y());
	}
	*item = new Item(coordinates);
}


void _get_relative_traversals(int i, int j, int m, int n, int pi, Tuple *path,
		list<vector<Tuple*>*>& relative_traversals) {

    // Reached the bottom of the matrix so we are left with
    // only option to move right
    if (i == m - 1) {
        for (int k = j; k < n; k++) {
			path[pi + k - j].set_x(i);
			path[pi + k - j].set_y(k);
		}

		vector<Tuple*> *relative_path = new vector<Tuple*>;
        for (int l = 0; l < pi + n - j; l++) {
			Tuple *tuple = new Tuple(path[l].get_x(), path[l].get_y());
			relative_path->push_back(tuple);
		}
		relative_traversals.push_back(relative_path);
        return;
    }

    // Reached the right corner of the matrix we are left with
    // only the downward movement.
    if (j == n - 1) {
        for (int k = i; k < m; k++) {
			path[pi + k - i].set_x(k);
			path[pi + k - i].set_y(j);
		}

		vector<Tuple*> *relative_path = new vector<Tuple*>;
        for (int l = 0; l < pi + m - i; l++) {
			Tuple *tuple = new Tuple(path[l].get_x(), path[l].get_y());
			relative_path->push_back(tuple);
		}
		relative_traversals.push_back(relative_path);
        return;
    }

    // Add the current cell to the path being generated
    //path[pi] = mat[i*n + j];
	path[pi].set_x(i);
	path[pi].set_y(j);

	int next_i = i + 1;
	int next_j = j + 1;

    //get all the paths that are possible after moving down
	if (next_i     == j*m/n || j ==  next_i*n/m ) {
		//|| next_i + 1 == j*m/n || j == (next_i + 1)*n/m
		//|| next_i - 1 == j*m/n || j == (next_i - 1)*n/m )  {
    	_get_relative_traversals(i+1, j, m, n, pi + 1, path, relative_traversals);
	}

    //get all the paths that are possible after moving right
	if (i     == next_j*m/n || next_j == i*n/m ) {
		//|| i + 1 == next_j*m/n || next_j == (i + 1)*n/m
		//|| i - 1 == next_j*m/n || next_j == (i - 1)*n/m) {
    	_get_relative_traversals(i, j+1, m, n, pi + 1, path, relative_traversals);
	}

    //get all the paths that are possible after moving diagonal
	if (next_i     == next_j*m/n || next_j == next_i*n/m ) {
		//|| next_i + 1 == next_j*m/n || next_j == (next_i + 1)*n/m
		//|| next_i - 1 == next_j*m/n || next_j == (next_i - 1)*n/m) {
    	_get_relative_traversals(i+1, j+1, m, n, pi + 1, path, relative_traversals);
	}
}

// The main function that gets all paths from top left to bottom right
// in a matrix of size mXn
void find_relevant_traversals(int m, int n, list<vector<Tuple*>*>& relative_traverals) {

	Tuple *path = new Tuple[m + n];
    _get_relative_traversals(0, 0, m, n, 0, path, relative_traverals);
	delete[] path;
}


void random_matrix(int K, int d, float **G, float from, float to) {
 	random_device rd{};
    mt19937 gen{rd()};
	normal_distribution<float> dis(0.5, 1);

	float U, V, S;
	for (size_t i = 0; i < K; i++) {
		for (size_t j = 0; j < d;) {
			double x = dis(gen);
			if (x > from && x < to) {
				G[i][j] = x;
				j++;
			}
		}
	}
}

void read_command_line_arguments(char *argv[], int& argc, string& input_file, string& query_file,
	string& output_file, int& k, int& L, int& w, int& st, float& eps, int& M_table) {
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
		else if (strcmp(argv[i], "-–k_vec") == 0) {
			k = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-L_grid") == 0) {
			L = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-w") == 0) {
			w = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-st") == 0) {
			st = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "--eps") == 0) {
			eps = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-M") == 0) {
			M_table = atoi(argv[i + 1]);
		}
	}
}

void read_command_line_arguments(char *argv[], int& argc, string& input_file, string& query_file,
	string& output_file, int& k, int& L, int& w, int& st, bool& check_for_identical_grid_flag, int& delta) {
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
		else if (strcmp(argv[i], "-–k_vec") == 0) {
			k = atoi(argv[i + 1]);
		}
		else if (strcmp(argv[i], "-L_grid") == 0) {
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
		else if (strcmp(argv[i], "--delta") == 0) {
			delta = atoi(argv[i + 1]);
		}
	}
}

void read_command_line_arguments(char *argv[], int& argc, string& input_file, string& query_file,
	string& output_file, int& k, int& L, int& w, int& st) {
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
	}
}

void delete_curves(list<Curve*> curves) {
	for (Curve *curve : curves) {
		delete curve;
	}
}

int read_2d_curves_from_file(string file_name, list<Curve*>& curves, int& max_length, int M_table) {
	string line, coordinate, name, tuple, part1, part2;
	int length;
	float x, y;
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
			if (length > M_table) {
				iss.ignore();
				continue;
			}

			points = new vector<Point*>;
			while ( iss >> part1 >> part2) {
				tuple = part1 + part2;
				sscanf(tuple.c_str(), "(%f, %f)", &x, &y);
				point = new Point();
				point->insert_coordinate(x);
				point->insert_coordinate(y);
				points->push_back(point);
			}
			curve = new Curve(name, points);
			curves.push_back(curve);
			//cout << "READING \n";
		}
		inputfile.close();
	}
	return 0;
}

int read_2d_curves_from_file(string file_name, list<Curve*>& curves, int& max_length) {
	string line, coordinate, name, tuple, part1, part2;
	int length;
	float x, y;
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
				point = new Point();
				point->insert_coordinate(x);
				point->insert_coordinate(y);
				points->push_back(point);
			}
			curve = new Curve(name, points);
			curves.push_back(curve);
			//cout << "READING \n";
		}
		inputfile.close();
	}
	return 0;
}

void print_curves(list<Curve*> curves)  {
	for (Curve *curve : curves) {
		curve->print();
		cout <<endl;
	}
}

void print_range_results(list<Item*> items, float radious)  {
	cout <<"R-near neighbors :"<<radious<<endl;
	for (Item *item : items) {
		cout << item->get_name();
		cout <<endl;
	}
}
void print_range_results_to_file(list<Item*> items,FILE *out ,float range){
	fprintf(out,"R-near neighbors : %lf", range);
	for (Item *item : items)
		fprintf(out,  "%s\n",   item->get_name().c_str() );
}

double DTW(vector<Point*>& p, vector<Point*>& q) {
	int m1 = p.size();
	int m2 = q.size();

	double dtw_array[m1][m2];

	//initialize (0,0) point of the array
	dtw_array[0][0] = euclidean_distance_2d(p[0], q[0]);

	//cout <<"p_x:"<<p[0]->get_x()<<endl;
	//cout <<"p_y:"<<p[0]->get_y()<<endl;
	//cout <<"q_x:"<<q[0]->get_x()<<endl;
	//cout <<"q_y:"<<q[0]->get_y()<<endl;
//
	//cout <<"diff:"<<p[0]->get_x() - q[0]->get_x()<<endl;
//
//
	//cout <<"dtw00: "<<dtw_array[0][0]<<endl;

	//initialize first row of the array
	for (size_t j = 1; j < m2; j++) {
		dtw_array[0][j] = euclidean_distance_2d( p[0], q[j] ) + dtw_array[0][j - 1];
	}

	//initialize first column of the array
	for (size_t i = 1; i < m1; i++) {
		dtw_array[i][0] = euclidean_distance_2d( p[i], q[0] ) + dtw_array[i - 1][0];
	}

	//calculate the rest cells of the array( (1,1) is calculated twice, but that's ok )
	for (size_t i = 1; i < m1; i++) {
		for (size_t j = 1; j < m2; j++) {
			dtw_array[i][j] = euclidean_distance_2d( p[i], q[j] ) +
				min({dtw_array[i - 1][j - 1], dtw_array[i][j - 1], dtw_array[i - 1][j]});
		}
	}
	return dtw_array[m1 - 1][m2 - 2];
}

float manhattan_distance_2d(Point *p, Point *q) {
	return abs(p->get_x() - q->get_x()) + abs(p->get_y() - q->get_y());
}

float euclidean_distance_2d(Point *p, Point *q) {
	return (p->get_x() - q->get_x())*(p->get_x() - q->get_x() ) +
	(p->get_y() - q->get_y())*(p->get_y() - q->get_y());
}

unsigned g_hash_function(vector<Type> x , int dimension, int w, int k,
	int bits_of_each_hash, unsigned M, vector<vector<float>*>& s_array, vector<unsigned>& m_powers) {

	unsigned hash_value, total_hash_value = 0;
	for (int i = 0; i < k; ++i) {
		hash_value = 0;
		hash_value = hash_function(x, dimension, w, M, *(s_array[i]), m_powers);
		//cout <<"hash: "<<hash_value<<endl;
		total_hash_value |= hash_value << (32 - (i + 1)*bits_of_each_hash);
	}

	//cout <<"total_hash_value: "<<total_hash_value<<endl;

	return total_hash_value;
}

unsigned hash_function(vector<Type> x, int dimension, int w, unsigned M,
	vector<float>& s, vector<unsigned>& m_powers) {

	vector<int> a;

	for (size_t i = 0; i < dimension; i++) {
		a.push_back( floor( (x[i] - s[i]) / w) );
	}

	unsigned sum = 0;
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

void print_results(Query_Result ann_result,string type,Query_Result exhaustive_result){
	if (ann_result.get_name() == "NULL") {
		cout <<"NOT FOUND"<<endl;
		return;
	}
	printf("Query: %s \nNearest neighbor: %s \ndistance%s: %f \ndistanceTrue: %f \n"
		"t%s: %f \n\n\n",ann_result.get_name().c_str(),exhaustive_result.get_name().c_str(),
		type.c_str(),ann_result.get_best_distance(),exhaustive_result.get_best_distance(),
		type.c_str(),ann_result.get_time());
}

void print_results(Query_Result ann_result,string type,string hashing,Query_Result exhaustive_result){

	printf("Query: %s \nMethod: %s \nHashFunction: %s\nFound Nearest: %s",ann_result.get_name().c_str(),type.c_str(),hashing.c_str());
}

void print_results_to_file(Query_Result ann_result,string type,FILE *out,Query_Result exhaustive_result){

	if (out==NULL){
		fprintf(stderr, "Error opening file '%s'\n", optarg);
	}
	//printf("Query: %s \nNearest neighbor: %s \ndistance%s: %ld\ndistanceTrue: %ld\nt%s: %ld \n\n\n",ann_result.get_name().c_str(),exhaustive_result.get_name().c_str(),type.c_str(),ann_result.get_best_distance(),exhaustive_result.get_best_distance(),type.c_str(),ann_result.get_time());
	fprintf(out, "Query: %s \nNearest neighbor: %s \ndistance%s: %f \ndistanceTrue: %f \n"
	"t%s: %f \n\n\n",ann_result.get_name().c_str(),exhaustive_result.get_name().c_str(),
	type.c_str(),ann_result.get_best_distance(),exhaustive_result.get_best_distance(),
	type.c_str(),ann_result.get_time());
	//fclose(out);
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

void print_parameters(int L, int k, int w, int search_threshold) {
	cout <<"L: "<<L<<endl;
	cout <<"k: "<<k<<endl;
	cout <<"w: "<<w<<endl;
	cout <<"search_threshold: "<<search_threshold<<endl;
}

void print_parameters(int L, int k, int w, int search_threshold, int dimension, int delta) {
	cout <<"delta: "<<delta<<endl;
	print_parameters(L, k, w, search_threshold, dimension);
}

void print_parameters(int L, int k, int w, int search_threshold, int dimension, float range) {
	cout <<"range: "<<range<<endl;
	print_parameters(L, k, w, search_threshold, dimension);
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

void read_vectors_from_file(string file_name, list<Item*>& items, float& radious) {
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
			if (name == "Radius:") {
				iss >> radious;
				continue;
			}
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
	double best_distance = numeric_limits<double>::max();
	string best = "";

	time_t time;
	time = clock();
	for (Curve *curve : *curves) {
		double cur_distance = DTW(query->get_points(), curve->get_points());
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

uint64_t  mul_mod2(uint64_t  a, uint64_t  b, uint64_t  m) {
	long double x;
	uint64_t c;
	int64_t  r;
	if (a >= m) a %= m;
	if (b >= m) b %= m;
	x = a;
	c = x * b / m;
	r = (int64_t)(a * b - c * m) % (int64_t)m;
	return r < 0 ? r + m : r;
}

uint64_t pow_mod(uint64_t  a, uint64_t  b, uint64_t  m) {
	uint64_t r = m==1?0:1;
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

int hammingDistance(unsigned n1, unsigned n2) {
    unsigned x = n1 ^ n2;
    int setBits = 0;

    while (x > 0) {
        setBits += x & 1;
        x >>= 1;
    }

    return setBits;
}
