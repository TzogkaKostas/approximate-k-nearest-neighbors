#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>
#include <time.h>
#include <cmath>
#include <limits>

using namespace std;

#include "../point/point.hpp"
#include "../curve/curve.hpp"
#include "../hash_table/hash_table.hpp"
#include "../helping_functions/helping_functions.hpp"
#include "../query_result/query_result.hpp"
#include "../Tuple/tuple.hpp"
#include "../lsh/lsh.hpp"
#include "../curve_projection_lsh/curve_projection_lsh.hpp"

#define L_DEFAULT 1
#define K_DEFAULT 4
#define W_DEFAULT 100
#define SEARCH_THRESHOLD (L_DEFAULT*10)
#define TABLE_SIZE_DIVIDED_BY 16
#define CURVE_DIMENSION_DEFAULT 2
#define CHECK_FOR_IDENTICAL_GRID_FLAG_DEFAULT false
#define EPS_DEFAULT 0.5


int main(int argc, char *argv[]) {
	/*
	int mat[5][5] = {
		{11, 12, 13, 14, 15},
		{21, 22, 23, 24, 25},
		{31, 32, 33, 34, 35},
		{41, 42, 43, 44, 45},
		{51, 52, 53, 54, 55} };
	int mat[3][3] = {
		{11, 12, 13},
		{21, 22, 23},
		{31, 32, 33}
	};
	int mat[2][2] = {
		{11, 12},
		{21, 22}
	};
	list<vector<Tuple*>*> rel_list;
    printAllPaths(*mat, 2, 2, rel_list);
	list<vector<Tuple*>*> relevant_traversals;
	get_relative_traversals(3, 3, relevant_traversals);
	cout << relevant_traversals.size();
	return 0;
	*/

	int L = L_DEFAULT;
	int k = K_DEFAULT;
	int w = W_DEFAULT;
	int search_threshold = SEARCH_THRESHOLD;
	int table_size_divived_by = TABLE_SIZE_DIVIDED_BY;
	int curve_dimension = CURVE_DIMENSION_DEFAULT;
    int eps = EPS_DEFAULT;

	//if (argc < 7 ) {
	//	cout <<"usage: ./lsh –d <input file> –q <query file> –k <int> -L <int> -ο <output file>"<<endl;
	//	return 1;
	//}

	//READ COMMAND LINE ARGUMENTS
	string input_file, query_file, output_file;
	read_command_line_arguments(argv, argc, input_file, query_file, output_file,
			k, L, w, search_threshold, eps, M_table);

	//READ CURVES FROM THE INPUT FILE
	list<Curve*> input_curves;
	int max_curve_length;
	read_2d_curves_from_file(input_file, input_curves, max_curve_length);

	//INITIALIZE PARAMETERS
	int table_size = max( (int)floor(input_curves.size()) /
		table_size_divived_by, table_size_divived_by );
	int hash_table_dimension = curve_dimension*max_curve_length; //<---------------------------------
	unsigned m = numeric_limits<unsigned>::max() + 1 - 5;
	search_threshold = max(10*L, search_threshold);
    int K_matrix = 0 - curve_dimension*log(eps/(eps*eps));

	cout <<"K_matrix: "<<K_matrix<<endl;
	cout <<"M: "<<M_table<<endl;
	print_parameters(L, k, w, search_threshold, hash_table_dimension);
	cout <<"K_matrix: "<<K_matrix<<endl;
	cout <<"Max curve length: "<<max_curve_length<<endl;

	//CREATE THE GRID STRUCTURE FOR CURVES WITH LSH
	Curve_Projection_LSH grid_projection(L, hash_table_dimension, w, k,
		curve_dimension, m, M_table, K_matrix);

	//INSERT INPUT DATA
	time_t time = clock();
	for(Curve *curve : input_curves) {
		grid_projection.insert_curve(curve);
	}
	time = clock() - time;
	cout <<"Data insertion time: "<< ((double)time) / CLOCKS_PER_SEC <<endl<<endl;

	delete_curves(input_curves);
	cout <<"main"<<input_curves.size()<<endl;
	return 0;

	//READ QUERY CURVES FROM THE INPUT FILE
	list<Curve*> queries;
	read_2d_curves_from_file(query_file, queries, max_curve_length);

	//HANDLE QUERIES
	Query_Result ann_query_result, exhaustive_query_result;
	time = clock();
	double sum_query_time = 0;
	double max_rate = -1;
	double sum_rate = 0;
	int found_nearest = 0;
	int total_distances = 0;
	int not_null = 0;
	for(Curve *query: queries) {
		cout <<"Query:"<<query->get_name()<<endl;

		//approximate nearest neighbor
		grid_projection.ANN(query, search_threshold, ann_query_result);
		print_ann_results(ann_query_result);

		//Exact nearest neighbor
		exhaustive_curve_search(&input_curves, query, exhaustive_query_result);
		print_exhaustive_search_results(exhaustive_query_result);
		cout <<"-------------------------------------------------------"<<endl;
		cout<<endl;

		//statistics info for searches that succeeded
		if (ann_query_result.get_time() != -1 && 0) {
			sum_query_time += ann_query_result.get_time();
			if (exhaustive_query_result.get_best_distance() != 0) { //division by zero
				max_rate = max(max_rate,
						(double)ann_query_result.get_best_distance()/exhaustive_query_result.get_best_distance());
				sum_rate += ann_query_result.get_best_distance()/exhaustive_query_result.get_best_distance();
			}
			if (ann_query_result.get_name() == exhaustive_query_result.get_name()) {
				found_nearest++;
			}
			total_distances += ann_query_result.get_best_distance();
			not_null++;
		}
	}
	time = clock() - time;
	cout <<endl;
	cout <<"Handling of queries time: "<< ((double)time) / CLOCKS_PER_SEC<<endl;
	cout << "Average query time: "<<sum_query_time/not_null<<endl;
	cout << "Max rate: "<<max_rate<<endl;
	cout << "Average rate: "<<sum_rate/not_null<<endl;
	cout << "Found "<<not_null<<"/"<<queries.size()<<" approximate nearest neighbors"<<endl;
	cout << "Found "<<found_nearest<<"/"<<queries.size()<<" exact nearest neighbors"<<endl;
	cout << "Average distance: "<<total_distances/queries.size()<<endl;

	delete_curves(input_curves);
	delete_curves(queries);
	return 0;
}
