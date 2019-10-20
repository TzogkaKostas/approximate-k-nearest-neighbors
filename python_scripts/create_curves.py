import random
import numpy as np
import sys

#GENERATE RANDOM INPUT CURVES
number_of_curves = int(sys.argv[1])
max_curve_length = 40
start_random = -500.0
end_random = 500.0

curves = []
for i in range(0, number_of_curves):
	points = []
	cur_max_length = np.random.uniform(1, max_curve_length, 1)
	for j in range(0, cur_max_length):
		point = np.random.uniform(start_random, end_random, 2)
		points.append(tuple(point))
	curves.append(points)

outfile = open('../grid_lsh_main/in.txt', 'w+')

i = 0
max_input_length = - 1
for c in curves:
	max_input_length = max(max_input_length, len(c))
	outfile.write("%d %d " %(i, len(c)))
	outfile.write(' '.join(map(str, c)))
	outfile.write("\n")
	#print ' '.join(map(str, c))
	i = i + 1


#GENERATE RANDOM QUERY CURVES
number_of_curves = int(sys.argv[2])
start_random = -500.0
end_random = 500.0

curves = []
for i in range(0, number_of_curves):
	points = []
	cur_max_length = np.random.uniform(1, max_input_length, 1)
	for j in range(0, cur_max_length):
		point = np.random.uniform(start_random, end_random, 2)
		points.append(tuple(point))
	curves.append(points)

outfile = open('../grid_lsh_main/q.txt', 'w+')
i = 0
max_input_length = - 1
for c in curves:
	max_input_length = max(max_input_length, len(c))
	outfile.write("%d %d " %(i, len(c)))
	outfile.write(' '.join(map(str, c)))
	outfile.write("\n")
	#print ' '.join(map(str, c))
	i = i + 1
