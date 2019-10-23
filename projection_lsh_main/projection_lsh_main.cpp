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
#define W_DEFAULT 400
#define SEARCH_THRESHOLD (L_DEFAULT*100)
#define CURVE_DIMENSION_DEFAULT 2
#define CHECK_FOR_IDENTICAL_GRID_FLAG_DEFAULT false
#define EPS_DEFAULT 0.5
#define M_TABLE_DEFAULT 4


int main(int argc, char *argv[]) {
	int L = L_DEFAULT;
	int k = K_DEFAULT;
	int w = W_DEFAULT;
	int search_threshold = SEARCH_THRESHOLD;
	int curve_dimension = CURVE_DIMENSION_DEFAULT;
    float eps = EPS_DEFAULT;
	int M_table = M_TABLE_DEFAULT;

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
	read_2d_curves_from_file(input_file, input_curves, max_curve_length, M_table);

	//INITIALIZE PARAMETERS
	unsigned m = numeric_limits<unsigned>::max() + 1 - 5;
	search_threshold = max(100*L, search_threshold);
    int K_matrix = 0 - curve_dimension*log2(eps)/(eps*eps);

	cout <<"K_matrix: "<<K_matrix<<endl;
	cout <<"M: "<<M_table<<endl;
	print_parameters(L, k, w, search_threshold);
	cout <<"Max curve length: "<<max_curve_length<<endl;
	cout <<"num of curves: "<<input_curves.size()<<endl;

	//CREATE THE GRID STRUCTURE FOR CURVES WITH LSH
	Curve_Projection_LSH grid_projection(L, w, k,
		curve_dimension, m, M_table, K_matrix);

	grid_projection.print_hash_tables();


	//INSERT INPUT DATA
	time_t time = clock();
	for(Curve *curve : input_curves) {
		grid_projection.insert_curve(curve);
	}
	time = clock() - time;
	cout <<"Data insertion time: "<< ((double)time) / CLOCKS_PER_SEC <<endl<<endl;

	//READ QUERY CURVES FROM THE INPUT FILE
	list<Curve*> queries;
	read_2d_curves_from_file(query_file, queries, max_curve_length, M_table);

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
		if (ann_query_result.get_time() != -1 ) {
			sum_query_time += ann_query_result.get_time();
			cout <<"kosdas"<<endl;
			if (exhaustive_query_result.get_best_distance() != 0) { //division by zero
				cout <<"kosdas"<<endl;
				max_rate = max(max_rate,
						(double)ann_query_result.get_best_distance()/exhaustive_query_result.get_best_distance());
				sum_rate += ann_query_result.get_best_distance()/exhaustive_query_result.get_best_distance();
			}
			cout << ann_query_result.get_name()<<endl;
			cout << exhaustive_query_result.get_name()<<endl;
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
