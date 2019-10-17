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
#include "../grid_projection/grid_projection.hpp"
#include "../hash_table/hash_table.hpp"
#include "../helping_functions/helping_functions.hpp"
#include "../query_result/query_result.hpp"
#include "../lsh/lsh.hpp"

#define L_DEFAULT 5
#define K_DEFAULT 4
#define W_DEFAULT 100
#define SEARCH_THRESHOLD (L_DEFAULT*3)
#define TABLE_SIZE_DIVIDED_BY 16 
#define CURVE_DIMENSION_DEFAULT 2
#define CHECK_FOR_IDENTICAL_GRID_FLAG_DEFAULT false


int main(int argc, char *argv[]) {
	int L = L_DEFAULT;
	int k = K_DEFAULT;
	int w = W_DEFAULT;
	int search_threshold = SEARCH_THRESHOLD;
	int table_size_divived_by = TABLE_SIZE_DIVIDED_BY;
	int curve_dimension = CURVE_DIMENSION_DEFAULT;
	bool check_for_identical_grid_flag = CHECK_FOR_IDENTICAL_GRID_FLAG_DEFAULT;

	//if (argc < 7 ) {
	//	cout <<"usage: ./lsh –d <input file> –q <query file> –k <int> -L <int> -ο <output file>"<<endl;
	//	return 1;
	//}

	//READ COMMAND LINE ARGUMENTS
	string input_file, query_file, output_file;
	read_command_line_arguments(argv, argc, input_file, query_file, output_file,
			k, L, w, search_threshold, check_for_identical_grid_flag);

	//READ CURVES FROM THE INPUT FILE
	list<Curve*> input_curves;
	int max_curve_length;
	read_2d_curves_from_file(input_file, input_curves, max_curve_length);

	//INITIALIZE PARAMETERS
	int table_size = max( (int)floor(input_curves.size()) / table_size_divived_by, table_size_divived_by );
	int hash_table_dimension = CURVE_DIMENSION_DEFAULT*max_curve_length;
	int delta = 4*CURVE_DIMENSION_DEFAULT*hash_table_dimension;
	unsigned m = numeric_limits<unsigned>::max() - 5;
	search_threshold = max(3*L, search_threshold);
	print_parameters(L, k, w, search_threshold, hash_table_dimension);

	//CREATE THE GRID STRUCTURE FOR CURVES WITH LSH 
	Grid_Projection grid_projection(L, hash_table_dimension, w, k, delta, curve_dimension, m);

	//INSERT INPUT DATA
	time_t time = clock();
	for(Curve *curve : input_curves) {
		grid_projection.insert_curve(curve);
	}
	time = clock() - time;
	cout <<"Data insertion time: "<< ((double)time) / CLOCKS_PER_SEC <<endl<<endl;
	

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
		grid_projection.ANN(query, search_threshold, ann_query_result, check_for_identical_grid_flag);
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
			if (ann_query_result.get_name() == exhaustive_query_result.get_name()) { //division by zero
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