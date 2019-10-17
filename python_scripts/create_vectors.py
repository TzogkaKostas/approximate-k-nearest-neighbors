import random
import numpy as np

number_of_vectors = 1
dimension = 4 + 1 # + 1 gia to name malaka ksupna
start_random = 0
end_random = 10

vectors = []
for i in range(0, number_of_vectors):
	coordinates = np.random.uniform(start_random, end_random, dimension)
	coordinates = map(int, coordinates)
	vectors.append(coordinates)

for v in vectors:
	print ' '.join(map(str, v))