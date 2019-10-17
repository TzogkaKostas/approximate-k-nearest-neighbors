import random
import numpy as np

number_of_curves = 5
max_curve_length = 600
start_random = -200.0
end_random = 200.0

curves = []
for i in range(0, number_of_curves):
	points = []
	for j in range(0, max_curve_length):
		point = np.random.uniform(start_random, end_random, 2)
		points.append(point)
	curves.append(points)

for v in curves:
	print(v)
	#print ' '.join(map(str, v))