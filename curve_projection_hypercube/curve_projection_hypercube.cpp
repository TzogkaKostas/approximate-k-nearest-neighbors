#include "curve_projection_hypercube.hpp"

using namespace std;

Curve_Projection_hypercube::Curve_Projection_hypercube(int hash_table_dimension, int w, int k, int delta,
        int curve_dimension, unsigned m,unsigned M ,int table_size,int probes){


}


Curve_Projection_hypercube::~Curve_Projection_hypercube() {
	for (size_t i = 0; i < L; i++) {
		for (size_t j = 0; j < L; j++) {
			delete table[i][j];
		}
	}
}
