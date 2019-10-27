// g++ hypercube.cpp hypercube_main.cpp ../helping_functions/helping_functions.cpp ../item/item_implem.cpp ../curve/curve_implem.cpp ../point/point_implem.cpp -g3
//./a.out -d trajectories_dataset -q query -o output
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <unordered_map>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdio>

using namespace std;
typedef float Type;

#include "../hypercube/hypercube.hpp"
#include "../helping_functions/helping_functions.hpp"
#include "../query_result/query_result.hpp"
#include "../item/item.hpp"
#define M_DEFAULT 500
#define K_DEFAULT 2
#define W_DEFAULT 4000
#define PROBES_DEFAULT 14


int main(int argc, char *argv[]) {
		int k;
		int k_s_g = K_DEFAULT;
		int w = W_DEFAULT;
		int M = M_DEFAULT;
		int probes = PROBES_DEFAULT;
		float radious = -1000;
		int PRINT_ON_SCREAN=0;

	    //READ COMMAND LINE ARGUMENTS
	  string input_file, query_file, output_file;
		int flag_defult=-1;
		read_command_line_arguments_hypercube(argv, argc, input_file, query_file,output_file,k,M,probes,flag_defult);

		if(output_file==""){
				PRINT_ON_SCREAN=1;
		}
	    //READ ITEMS FROM THE INPUT FILE
	  list<Item*> input_items;
	  read_vectors_from_file(input_file, input_items);
		int table_size=k;//initilized after insert items
		//w = calculate_w(input_items);

	  int dimension = input_items.front()->get_coordinates()->size();
		//	table_size=log(input_items.size());
		if (flag_defult==-1){
			table_size=log2(input_items.size());
		}
    unsigned m = numeric_limits<unsigned>::max() + 1 - 5;
    //print_parameters(0, table_size, w, probes, dimension);
		cout << "k " << k_s_g<<endl;
		cout << "M " << M<<endl;
		cout << "probes " << probes << endl;
		//CREATE THE HYPERCUBE STRUCTURE

		Hypercube hypercube (table_size,dimension,w,k_s_g,m,M);
		//cout << "ARGS " << input_file << " " << query_file << " " << output_file ;

		//INSERT INPUT DATA
		time_t time = clock();
		for(Item *item: input_items) {
			hypercube.insert_item(item);
		}
		// unordered_multimap<unsigned, Item*>* print;
		// print =hypercube.hash_table->get_f_values_map();
		// unordered_multimap<unsigned, Item*>::iterator it = print->begin();
		// for (; it != print->end(); it++)
	  //       cout << "<" << it->first << ", " << it->second->get_name()
	  //            << ">  ";
	  // unsigned nbuckets = hypercube.hash_table->get_f_values_map()->bucket_count();
		// std::cout << "myumm has " << nbuckets << " buckets:\n";
		//
	  // for (unsigned i=0; i<nbuckets; i++) {
		//   for (auto it = hypercube.hash_table->get_f_values_map()->begin(i);it!= hypercube.hash_table->get_f_values_map()->end(i);it++){
	  //     	std::cout << "[" << hypercube.hash_table->get_f_values_map()->begin(i)->first << ":" << it->second->get_name() << "] ";
		//   	break;
	  // 		}
	  // }
		// cout << endl;

		time = clock() - time;
		cout <<"Data insertion total time: "<< ((double)time) / CLOCKS_PER_SEC <<endl<<endl;

		//HANDLE QUERIES
		list<Item*> queries;
		read_vectors_from_file(query_file, queries, radious);
		cout <<"Radious: "<<radious<<endl;



		Query_Result ann_query_result, exhaustive_query_result,range_query_result;
		time = clock();
		double sum_query_time = 0;
		double max_rate = -1;
		double sum_rate = 0;
		int found_nearest = 0;
		int total_distances = 0;
		int not_null = 0;
		list<Item*> range_items;
		FILE *out;
		out= fopen(output_file.c_str(), "w");
		for(Item *query: queries) {
			//cout <<"Query:"<<query->get_name()<<endl;
			//approximate nearest neighbor
			hypercube.ANN(query, probes, ann_query_result);
			//print_ann_results(ann_query_result);

			//Exact nearest neighbor
			exhaustive_search(&input_items, query, exhaustive_query_result);
			//print_exhaustive_search_results(exhaustive_query_result);
			if(PRINT_ON_SCREAN==1)
					print_results(query->get_name(),ann_query_result,"Cube", exhaustive_query_result);
			else
					print_results_to_file(query->get_name(),ann_query_result,"Cube",out ,exhaustive_query_result);

			//cout<<endl;
			if (radious > 0) {
				hypercube.range_search(query, probes, radious, range_items, range_query_result);
				print_range_results(range_items, radious);
				range_items.clear();
			}

			//cout <<"--------------------------------------------------------"<<endl;
			//statistics info
			if (ann_query_result.get_time() != -1) {
					sum_query_time += ann_query_result.get_time();
					max_rate = max(max_rate, (double)ann_query_result.get_best_distance()/exhaustive_query_result.get_best_distance());
					sum_rate += ann_query_result.get_best_distance()/exhaustive_query_result.get_best_distance();
					if (ann_query_result.get_name() == exhaustive_query_result.get_name()) {
							found_nearest++;
					}
					total_distances += ann_query_result.get_best_distance();
					not_null++;
			}
			//cout <<"\nEND \n";

		}
		time = clock() - time;
		cout <<"Handling of queries(ann and enn) total time: "<< ((double)time) / CLOCKS_PER_SEC<<endl;
		cout << "Average query time: "<<sum_query_time/not_null<<endl;
		cout << "Max rate: "<<max_rate<<endl;
		cout << "Average rate: "<<sum_rate/not_null<<endl;
		cout << "Found "<<not_null<<"/"<<queries.size()<<" approximate nearest neighbors"<<endl;
		cout << "Found "<<found_nearest<<"/"<<queries.size()<<" exact nearest neighbors"<<endl;
		cout << "Average distance: "<<total_distances/queries.size()<<endl;
		cout <<"--------------------------------------------------------------------------"<<endl;
		cout <<endl;

		delete_items(input_items);
		delete_items(queries);
		return 0;
}
