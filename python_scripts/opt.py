from subprocess import call

input_file = "../small_data/input_small.txt"
query_file = "../small_data/query_small.txt"

#best w: 5000, max=3.2, avg = 1.05
#best w: 6400, max=4.7, avg = 1.09
#best w: 6600, max=4.2, avg = 1.09


for w in range(4900, 5000, 100):
	for L in [5]:
		for k in [4]:
			for div in [8]:		
				for st in [1000]:
					call(["./main", "-d", input_file, "-q", query_file, "-L",
						str(L), "-k", str(k), "-st", str(st), "-div", str(div), "-w", str(w)])
