import subprocess
import sys
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns

programs = ['matrix.out', 'matrix.py']

names_progs = {'matrix.out':'c++', 'matrix_opt.out':'c++ optimized', 'matrix.py':'python3', 'matrix_opt.out':'c++ optimized'}

functions = ['block', 'ijk', 'ikj']

#sizes = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]

sizes = [64, 128, 256]

times = {}

ms_time_re = re.compile(r"\d*(\.)?\d*ms")

def plot_graph(times):

	sns.set(style='darkgrid')

	for program in programs:
		example_prg = times[program]
		for function in functions:
			lin_dic = {'size':example_prg[function][0], 'time':example_prg[function][1]}
			ax = sns.lineplot(x='size', y='time', legend='full', label=names_progs[program]+" "+function, markers=True, marker="o", data=lin_dic)
			ax.set(xlabel='size', ylabel='time')
		
		name = "times_for_imp_" + program + ".png"
		plt.legend()
		plt.savefig(name)

	plt.legend()
	plt.savefig("times_all_imps.png")
	print(times)
	plt.show()


for program in programs:

	program_loc = os.getcwd()+'/'+program

	for function in functions:
		for size in sizes:
			run_cmd = [program_loc, function, str(size)]
			
			if(program.split('.')[-1] == 'py'):
				run_cmd = ['python3'] + run_cmd

			result = subprocess.run(run_cmd, stdout=subprocess.PIPE)
			result_out = result.stdout.decode('utf-8')

			ms_time = float(ms_time_re.search(result_out)[0].split('ms')[0])

			if(times.get(program)==None):
				times[program] = {function:[[size], [ms_time]]}
			else:
				if(times.get(program).get(function) != None):
					times[program][function][0].append(size)
					times[program][function][1].append(ms_time)
				else:
					times[program][function] = [[size], [ms_time]]

			print('Running ' + function + " implementation in " + names_progs[program] + '.\n' + result_out)

plot_graph(times)

#print(times)