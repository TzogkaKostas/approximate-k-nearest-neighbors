//g++ grid_hypercube_main.cpp ../curve_grid_hypercube/curve_grid_hypercube.cpp ../item/item_implem.cpp ../curve/curve_implem.cpp ../query_result/query_result.hpp ../point/point_implem.cpp ../hash_table/hash_table_implem.cpp ../helping_functions/helping_functions.cpp ../relevant_traversals/relevant_traversals_implem.cpp ../Tuple/tuple.hpp -g3
//./a.out -d data -q query -o output
#include "../curve_grid_hypercube/curve_grid_hypercube.hpp"
#include "../item/item.hpp"
#include "../curve/curve.hpp"
#include "../query_result/query_result.hpp"
#include "../point/point.hpp"
#include "../hash_table/hash_table.hpp"
#include "../helping_functions/helping_functions.hpp"
#include "../relevant_traversals/relevant_traversals.hpp"
#include "../Tuple/tuple.hpp"

using namespace std;
#define M_DEFAULT 500
#define K_DEFAULT 4
#define W_DEFAULT 50
#define PROBES_DEFAULT 14
#define L_DEFAULT 4
#define CHECK_FOR_IDENTICAL_GRID_FLAG_DEFAULT false
#define CURVE_DIMENSION_DEFAULT 2
//k =2

int main(int argc, char *argv[]) {
    int k;
    int k_s_g = K_DEFAULT;
    int w = W_DEFAULT;
    int M = M_DEFAULT;
    int probes = PROBES_DEFAULT;
    int L=L_DEFAULT;
    float delta = -1;
    bool check_for_identical_grid_flag = CHECK_FOR_IDENTICAL_GRID_FLAG_DEFAULT;
    int curve_dimension = CURVE_DIMENSION_DEFAULT;
    int PRINT_ON_SCREAN=0;

    //READ COMMAND LINE ARGUMENTS
    string input_file, query_file, output_file;
    int flag_defult=-1;
    read_command_line_arguments_hypercube_grid(argv, argc, input_file, query_file,output_file,k,M,probes,L,delta,flag_defult);

    if(output_file==""){
				PRINT_ON_SCREAN=1;
		}
    //READ ITEMS FROM THE INPUT FILE
    list<Curve*> input_curves;
    int max_curve_length;
    int table_size=k;//initilized after insert items
    read_2d_curves_from_file(input_file, input_curves, max_curve_length);
    //cout << "input_curves : "<<input_curves.size()<<endl;
    int hash_table_dimension = curve_dimension*max_curve_length;
    if (flag_defult==-1){
        table_size=log2(input_curves.size());
    }
    if(flag_defult==-1||flag_defult==0 || flag_defult ==1){
      //cout <<"aaaaaaaa\n";
      delta =  calculate_delta(input_curves);
    }
    unsigned m = numeric_limits<unsigned>::max() + 1 - 5;

    cout << "L "<< L<<endl;
    cout << "k "<< k_s_g<<endl;
    cout << "Probe "<< probes<<endl;
    cout << "delta " << delta<<endl;

    //CREATE THE HYPERCUBE STRUCTURE
    Curve_Grid_hypercube h_curve_grid(L,hash_table_dimension,w,k_s_g,delta,curve_dimension,m,M,table_size,probes);
    //cout << "ARGS " << input_file << " " << query_file << " " << output_file ;

    //INSERT INPUT DATA
    list<Curve*> grid_curves;
    time_t time = clock();
    for(Curve *curve : input_curves) {
        //cout <<curve->get_name()<<endl;
		h_curve_grid.insert_curve(curve, &grid_curves);
	}

    time = clock() - time;
	cout <<"Data insertion time: "<< ((double)time) / CLOCKS_PER_SEC <<endl<<endl;

    //grid_projection.print_hash_tables_names();

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
    FILE *out;
    out= fopen(output_file.c_str(), "w");
	for(Curve *query: queries) {

		//approximate nearest neighbor
		h_curve_grid.ANN(query, probes, ann_query_result, check_for_identical_grid_flag);

		//Exact nearest neighbor
		exhaustive_curve_search(&input_curves, query, exhaustive_query_result);

    if(PRINT_ON_SCREAN==1)
        print_results(query->get_name(),ann_query_result,"Grid","Cube", exhaustive_query_result);
    else
        print_results_to_file(query->get_name(),ann_query_result,"Grid","Cube",out ,exhaustive_query_result);

		//cout <<"-------------------------------------------------------"<<endl;
		//cout<<endl;

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
	cout << "Max rate: "<<max_rate<<endl;
	cout << "Average rate: "<<sum_rate/not_null<<endl;
	cout << "Found "<<not_null<<"/"<<queries.size()<<" approximate nearest neighbors"<<endl;
	cout << "Found "<<found_nearest<<"/"<<queries.size()<<" exact nearest neighbors"<<endl;
	cout << "Average distance: "<<total_distances/queries.size()<<endl;

	delete_curves(grid_curves);
	delete_curves(input_curves);
	delete_curves(queries);
	return 0;

}
