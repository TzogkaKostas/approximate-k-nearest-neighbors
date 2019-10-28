// g++ projection_hypercube_main.cpp ../item/item_implem.cpp ../curve/curve_implem.cpp ../point/point_implem.cpp ../helping_functions/helping_functions.cpp ../relevant_traversals_hypercube/relevant_traversals_hypercube.cpp ../curve_projection_hypercube/curve_projection_hypercube.cpp ../curve_grid_hypercube/curve_grid_hypercube.cpp -g3
#include <unordered_map>
#include <math.h>
#include <limits>
#include <random>
#include <iostream>
#include "../item/item.hpp"
#include "../curve/curve.hpp"
#include "../query_result/query_result.hpp"
#include "../point/point.hpp"
#include "../helping_functions/helping_functions.hpp"
#include "../relevant_traversals_hypercube/relevant_traversals_hypercube.hpp"
#include "../curve_projection_hypercube/curve_projection_hypercube.hpp"
#include "../Tuple/tuple.hpp"

#define M_DEFAULT 500
#define K_DEFAULT 4
#define W_DEFAULT 4000
#define PROBES_DEFAULT 2
#define CHECK_FOR_IDENTICAL_GRID_FLAG_DEFAULT false
#define CURVE_DIMENSION_DEFAULT 2
#define EPS_DEFAULT 0.5
#define M_TABLE_DEFAULT 4

int main(int argc, char *argv[]) {
    int k_s_g = K_DEFAULT;
    int table_size_hypercube;
	int w = W_DEFAULT;
    int M = M_DEFAULT;
    int probes = PROBES_DEFAULT;
    int curve_dimension = CURVE_DIMENSION_DEFAULT;
    float eps =EPS_DEFAULT;
    int M_table = M_TABLE_DEFAULT;
    int PRINT_ON_SCREAN=0;


    //READ COMMAND LINE ARGUMENTS
    string input_file, query_file, output_file;
    int flag_defult=-1;

    read_command_line_arguments_hypercube_projection(argv, argc, input_file, query_file,table_size_hypercube,M,probes,eps,output_file,flag_defult);
    if(output_file==""){
        PRINT_ON_SCREAN=1;
    }
    //READ ITEMS FROM THE INPUT FILE
    list<Curve*> input_curves;
    int max_curve_length;
    //int table_size=k;//initilized after insert items
    read_2d_curves_from_file(input_file, input_curves, max_curve_length,M_table);
    //cout << "input_curves : "<<input_curves.size()<<endl;
    int hash_table_dimension = curve_dimension*max_curve_length;
    if (flag_defult==-1){
        table_size_hypercube=log2(input_curves.size());
    }
    unsigned m = numeric_limits<unsigned>::max() + 1 - 5;
    print_parameters(1, table_size_hypercube, w, probes, hash_table_dimension);
    int K_matrix = 0 - curve_dimension*log2(eps)/(eps*eps);

    //CREATE THE GRID STRUCTURE FOR CURVES WITH HYPERCUBE
    Curve_Projection_hypercube grid_projection(w,k_s_g,curve_dimension,m,M,K_matrix,table_size_hypercube,probes,M_table);

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
    double total_distances = 0;
    int not_null = 0;
    FILE *out = fopen(output_file.c_str(), "w");
    for(Curve *query: queries) {
        //Approximate nearest neighbor
        grid_projection.ANN(query, probes, ann_query_result);

        //Exact nearest neighbor
        exhaustive_curve_search(&input_curves, query, exhaustive_query_result);

        if(PRINT_ON_SCREAN==1)
            print_results(query->get_name(),ann_query_result,"Projection","Cube", exhaustive_query_result);
        else {
            print_results_to_file(query->get_name(),ann_query_result,"Projection","Cube",out ,exhaustive_query_result);

        }

        //statistics info for searches that succeeded
        if (ann_query_result.get_time() != -1 ) {
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
    return 0;
}
