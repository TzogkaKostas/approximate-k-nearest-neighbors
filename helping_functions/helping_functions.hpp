#ifndef HELP_H
#define HELP_H

#include <list>
#include <vector>
#include <unordered_map>
#include <unistd.h>
#include <random>
#include "../query_result/query_result.hpp"
#include "../curve/curve.hpp"
#include "../point/point.hpp"
#include "../Tuple/tuple.hpp"

using namespace std;

unsigned mod(long long a, long long b);
unsigned add_mod(long long a, long long b, long long m);
uint32_t mul_mod(long long a, long long b, long long m);
uint64_t mul_mod2(uint64_t a, uint64_t b, uint64_t m) ;
uint64_t pow_mod(uint64_t a, uint64_t b, uint64_t m);
void random_float_list(double from, double to, list<float>& random_list, int size);
void random_float_vector(double from, double to, vector<float>& random_vector, int size);
int max_from_list(list<int> my_list);
int max_from_list(list<float> my_list);
int max_from_vector(vector<float> my_vector);
int max_from_vector(vector<int> my_vector);
int find_dimension_from_file(string file_name);
void print_vector(vector<int> my_list);
void print_vector(vector<float> my_list);
unsigned g_hash_function(vector<Type> x , int dimension, int w, int k,
	int bits_of_each_hash, unsigned M, vector<vector<float>*>& s_array, vector<unsigned>& m_powers);
unsigned hash_function(vector<Type> x, int dimension, int w, unsigned M,
	vector<float>& s, vector<unsigned>& m_powers);
void read_command_line_arguments(char *argv[], int& argc, string& input_file, string& query_file,
	string& output_file, int& k, int& L, int& w, int& st);
void read_command_line_arguments(char *argv[], int& argc, string& input_file, string& query_file,
	string& output_file, int& k, int& L, int& w, int& st, bool& check_for_identical_grid_flag, float& delta);
void read_command_line_arguments(char *argv[], int& argc, string& input_file, string& query_file,
	string& output_file, int& k, int& L, int& w, int& st, float& eps, int& M_table);
unsigned long long int manhattan_distance(vector<Type> x1, vector<Type> x2);
void print_ann_results(Query_Result ann_result);
void print_exhaustive_search_results(Query_Result exhaustive_result);
void print_parameters(int L, int k, int w, int search_threshold, int dimension);
void print_parameters(int L, int k, int w, int search_threshold);
void print_parameters(int L, int k, int w, int search_threshold, int dimension, float range);
void get_vector_from_line(string line, Item& item);
void read_vectors_from_file(string file_name, list<Item*>& items);
void read_vectors_from_file(string file_name, list<Item*>& items, float& range);
unsigned long long int manhattan_distance(vector<Type> x1, vector<Type> x2);
void exhaustive_search(list<Item*> *items, Item *query, Query_Result& query_result);
void delete_items(list<Item*> items);
int read_2d_curves_from_file(string file_name, list<Curve*>& curves, int& max_length, double& max_coord);
int read_2d_curves_from_file(string file_name, list<Curve*>& curves, int& max_length, int M_table);
void print_curves(list<Curve*> curves);
void delete_curves(list<Curve*> curves);
void exhaustive_curve_search(list<Curve*> *curves, Curve *query, Query_Result& query_result);
double DTW(vector<Point*>& p, vector<Point*>& q);
double manhattan_distance_2d(Point *p, Point *q);
double euclidean_distance_2d(Point *p, Point *q);
void snap_curve(Curve *curve, Point *t, Curve **snapped_curve, float delta);
void get_snapped_point(Point *point, float delta, Point *t, Point **snapped_point);
void zip_points(Curve *snapped_curve, Item **item);
void fill_curve(Curve *curve, int pad_length);
void convert_2d_curve_to_vector(Curve *curve, Point *t, double delta, int hash_table_dimesion,
	int curve_dimension, Curve **snapped_curve, Item **item, double max_coord);
void random_matrix(int K, int d, double **G, float mean, float deviation);
void print_range_results(list<Item*> items, float range);
void find_relevant_traversals(int m, int n, list<vector<Tuple*>*>& relative_traverals);
void print_range_results_to_file(list<Item*> items,FILE *out ,float range);
void convert_2d_curve_to_vector_by_projection(vector<Tuple*>& traversal, int U_or_V, double **G_matrix, Curve *curve,
		int G_rows, int G_cols, Item **item);
void matrix_multiplication(vector<Tuple*>& traversal, int U_or_V, double **G_matrix, Curve *curve,
		int G_rows, int G_cols, Item **item);
void read_command_line_arguments_hypercube(char *argv[], int& argc,string& input_file, string& query_file,
	string& output_file, int& k, int& M,int& probes,int &flag);
unsigned f_hash_function(vector<Type> x , int dimension,int w, int k,
		int bits_of_each_hash, unsigned M, vector<unsigned>& m_powers,vector<vector<float>*>& s_array,vector< unordered_map<unsigned,int> *>&g_value,int f);
int hammingDistance(unsigned n1, unsigned n2);
void read_command_line_arguments_hypercube_grid(char *argv[], int& argc,string& input_file, string& query_file,
	string& output_file, int& k, int& M,int& probes,int &L,double &delta ,int &flag );
void read_command_line_arguments_hypercube_projection(char* argv[],int &argc,string &input_file,string &query_file,int& k, int& M,int & probes,double &e,string &output_file,int &flag);
double calculate_delta(list<Curve*> curves);
int calculate_w(list<Item*> items);
void print_results(string query,Query_Result ann_result,string type,Query_Result exhaustive_result);
void print_results_to_file(string query,Query_Result ann_result,string type,FILE *out,Query_Result exhaustive_result);
void print_results(string query,Query_Result ann_result,string type,string hashing,Query_Result exhaustive_result);
void print_results_to_file(string query,Query_Result ann_result,string type,string hashing,FILE *out,Query_Result exhaustive_result);
void print_vector(vector<double> my_list);
void zip_points2(Curve *grid_curve, Item **item, int pad_length, double max_coord);
int calculate_curve_w(list<Curve*> curves);
#endif
