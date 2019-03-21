import sys
import pprint
import time
from copy import deepcopy


def ijk(A, B, C):
	for i in range(len(A)):
		for j in range(len(A)):
			for k in range(len(A)):
				C[i][j] += A[i][k] * B[k][j]


def ikj(A, B, C):
	for i in range(len(A)):
		for k in range(len(A)):
			for j in range(len(A)):
				C[i][j] += A[i][k] * B[k][j]


def block(A, B, C):
	return

def gen_matrix(size, gradual=True):
	matrix = []
	
	for i in range(size):
		line = []
		for j in range(size):
			line.append(i*size + j)
		matrix.append(line)

	return matrix

if __name__ == "__main__":

	if(len(sys.argv) < 3):
		print("Incorrect Number of arguments, type function and size of matrix.")

	else:

		size = int(sys.argv[2])
		imp = sys.argv[1]	

		A = gen_matrix(size)
		B = gen_matrix(size)
		C = gen_matrix(size)

		default_c = C.copy()


		t = time.process_time()
		elapsed_time = 0

		for i in range(10):
			locals()[imp](A, B, C)
			C = default_c.copy()

		elapsed_time = (time.process_time() - t)/10
		
		if( "-v" in sys.argv or "--verbose" in sys.argv):
			#print("Matrix A and B: {} \n Matrix C: {}".format(A, C))
			pprint.pprint("Matrix A and B: ")
			pprint.pprint(A)
			pprint.pprint("Matrix C: ")
			pprint.pprint(C)

		pprint.pprint("Time to run implementation: {0:0.3f}ms.".format(elapsed_time*1000))