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
#include "../curve_grid_lsh/curve_grid_lsh.hpp"
#include "../hash_table/hash_table.hpp"
#include "../helping_functions/helping_functions.hpp"
#include "../query_result/query_result.hpp"
#include "../lsh/lsh.hpp"

#define L_DEFAULT 5
#define K_DEFAULT 4
#define W_DEFAULT 40
#define SEARCH_THRESHOLD (L_DEFAULT*100)
#define CURVE_DIMENSION_DEFAULT 2
#define CHECK_FOR_IDENTICAL_GRID_FLAG_DEFAULT false


int main(int argc, char *argv[]) {
	int L = L_DEFAULT;
	int k = K_DEFAULT;
	int w = W_DEFAULT;
	int search_threshold = SEARCH_THRESHOLD;
	int curve_dimension = CURVE_DIMENSION_DEFAULT;
	bool check_for_identical_grid_flag = CHECK_FOR_IDENTICAL_GRID_FLAG_DEFAULT;
	float delta = -1;
	double max_coord;

	if (argc < 5 ) {
		cout <<"usage: ./lsh –d <input file> –q <query file> -ο <output file>"<<endl;
		return 1;
	}

	//READ COMMAND LINE ARGUMENTS
	string input_file, query_file, output_file;
	read_command_line_arguments(argv, argc, input_file, query_file, output_file,
			k, L, w, search_threshold, check_for_identical_grid_flag, delta);

	//READ CURVES FROM THE INPUT FILE
	list<Curve*> input_curves;
	int max_curve_length;
	read_2d_curves_from_file(input_file, input_curves, max_curve_length, max_coord);

	//INITIALIZE PARAMETERS
	int hash_table_dimension = curve_dimension*max_curve_length;
	if (delta == - 1) {
		delta = calculate_delta(input_curves);
	}
	if (w == -1) {
		w = calculate_curve_w(input_curves);
	}
	
	unsigned m = numeric_limits<unsigned>::max() + 1 - 5;
	search_threshold = max((int)input_curves.size()/10, search_threshold);
	print_parameters(L, k, w, search_threshold, hash_table_dimension);
	cout <<"delta: "<<delta<<endl;

	//CREATE THE GRID STRUCTURE FOR CURVES WITH LSH 
	Curve_Grid_LSH grid_projection(L, hash_table_dimension, w, k, delta, curve_dimension,
		m, max_coord);

	//INSERT INPUT DATA
	list<Curve*> grid_curves;
	time_t time = clock();
	for(Curve *curve : input_curves) {
		grid_projection.insert_curve(curve, &grid_curves);
	}
	time = clock() - time;
	cout <<"Data insertion time: "<< ((double)time) / CLOCKS_PER_SEC <<endl<<endl;

	//grid_projection.print_hash_tables_names();

	//READ QUERY CURVES FROM THE INPUT FILE
	list<Curve*> queries;
	read_2d_curves_from_file(query_file, queries, max_curve_length, max_coord);

	//HANDLE QUERIES
	Query_Result ann_query_result, exhaustive_query_result;
	time = clock();
	double sum_query_time = 0;
	double max_rate = -1;
	double sum_rate = 0;
	int found_nearest = 0;
	float total_distances = 0;
	int not_null = 0;
	FILE *out = fopen(output_file.c_str(), "w");
	for(Curve *query: queries) {
		//approximate nearest neighbor
		grid_projection.ANN(query, search_threshold, ann_query_result, check_for_identical_grid_flag);

		//Exact nearest neighbor
		exhaustive_curve_search(&input_curves, query, exhaustive_query_result);

		if (output_file != "") {
			print_results_to_file(query->get_name(), ann_query_result, "Grid", "LSH", out ,exhaustive_query_result);
		}
		else {
			print_results(query->get_name(), ann_query_result, "Grid", "LSH", exhaustive_query_result);
		}

		//statistics info for searches that succeeded
		if (ann_query_result.get_time() != -1) {
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
	cout << "Max AF: "<<max_rate<<endl;
	cout << "Average AF: "<<sum_rate/not_null<<endl;
	cout << "Found "<<not_null<<"/"<<queries.size()<<" approximate nearest neighbors"<<endl;
	cout << "Found "<<found_nearest<<"/"<<queries.size()<<" exact nearest neighbors"<<endl;
	cout << "Average distance: "<<total_distances/queries.size()<<endl;
	
	delete_curves(grid_curves);
	delete_curves(input_curves);
	delete_curves(queries);
	return 0;
}