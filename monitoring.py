import subprocess
import sys
import os
import re

programs = ['matrix.out', 'matrix.py']

names_progs = {'matrix.out':'c++', 'matrix.py':'python3', 'matrix_opt.out':'c++ optimized'}

functions = ['block', 'ijk', 'ikj']

#sizes = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]

sizes = [32, 64]

times = {}

ms_time_re = re.compile(r"\d*(\.)?\d*ms")

for program in programs:

	program_loc = os.getcwd()+'/'+program

	for function in functions:
		for size in sizes:
			run_cmd = [program_loc, function, str(size)]
			
			if(program.split('.')[-1] == 'py'):
				run_cmd = ['python3'] + run_cmd

			result = subprocess.run(run_cmd, stdout=subprocess.PIPE)
			result_out = result.stdout.decode('utf-8')

			

			ms_time = ms_time_re.search(result_out)[0]

			times[program] = {function:{size:ms_time}}

			print('Running ' + function + " implementation in " + names_progs[program] + '.\n' + result_out)

print(times)