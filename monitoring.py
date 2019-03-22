import subprocess
import sys
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
import copy

programs = ['matrix.out', 'matrix.py', 'SingleCore.class']

names_progs = {'matrix.out':'c++', 'SingleCore.class':'java', 'matrix_opt.out':'c++ optimized', 'matrix.py':'python3', 'matrix_opt.out':'c++ optimized'}

functions = ['block', 'ijk', 'ikj']

valgrind_call = "valgrind --tool=cachegrind " # valgrind --tool=cachegrind ./matrix.out

size1 = [x for x in range(600, 3001, 400)]
size2 = [x for x in range(4000, 10001, 2000)]

#sizes = size1 + size2

#sizes = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]

sizes = [32, 64, 128]

times = {}

ms_time_re = re.compile(r"\d*(\.)?\d*ms")

def plot_graph(times):

    sns.set(style='darkgrid')

    plots_backup = []

    for program in programs:
        example_prg = times[program]
        for function in functions:

            if(function == 'block' and program == 'matrix.py'):
                continue

            lin_dic = {'size':example_prg[function][0], 'time':example_prg[function][1]}
            plots_backup.append([lin_dic.copy(), function, program])
            ax = sns.lineplot(x='size', y='time', legend='full', label=names_progs[program]+" "+function, markers=True, marker="o", data=lin_dic)
            ax.set(xlabel='size', ylabel='time')
        
        name = "times_for_imp_" + program + ".png"
        plt.legend()
        plt.savefig(name)
        plt.clf()

    plt.legend()
    plt.clf()

    for plot in plots_backup:
        lin_dic, function, program = plot
        ax = sns.lineplot(x='size', y='time', legend='full', label=names_progs[program]+" "+function, markers=True, marker="o", data=lin_dic)
        ax.set(xlabel='size', ylabel='time')

    plt.savefig("times_all_imps.png")
    #plt.show()


def all_times():
    for program in programs:

        program_loc = os.getcwd()+'/'+program

        for function in functions:
            for size in sizes:
                run_cmd = [program_loc, function, str(size)]
                
                if(program.split('.')[-1] == 'py'):
                    run_cmd = ['python3'] + run_cmd

                if(program == 'SingleCore.class'):
                    run_cmd = ['java'] + run_cmd
                    print(run_cmd)


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


def cache_line():
    lscpu = subprocess.run(['lscpu'], stdout=subprocess.PIPE)
    lscpu_out = lscpu.stdout.decode('utf-8')

    l1d_re = re.compile(r"L1d\s*cache:\s*(\d*)K")
    l1d_size = int(l1d_re.search(lscpu_out)[1])

    l2_re = re.compile(r"L2\s*cache:\s*(\d*)K")
    l2_size = int(l2_re.search(lscpu_out)[1])

    print("L1d cache is: {}K of size.".format(l1d_size))
    print("L2 cache is:  {}K of size.".format(l2_size))

    return [l1d_size, l2_size]

def set_envsv_cache():

    with open("lscpu.txt", "w+") as output:
        subprocess.call(["lscpu"], stdout=output);

    L1, L2 = cache_line()

    print(str(L1))

    os.environ["L1DSIZE"] = str(L1)
    os.environ["L2SIZE"] = str(L2)

def call_valgrind(command):

    result = subprocess.run(run_cmd, stdout=subprocess.PIPE)
    result_out = result.stdout.decode('utf-8')

def cache_sizes_block():

    times = {}
    program = 'cpp'
    vals = [128, 256, 512]

    L1, L2 = cache_line()

    double_size = 8

    vals += [int(L1*1024//double_size)]
    vals += [int((L1*1024//double_size)**0.5)]

    plt.clf()

    for val in vals:
        for size in [128, 256]:
            env_var = 'L1DSIZE='+str(val)
            run_cmd = ['./matrix.out', 'block', str(size)]
            
            d = dict(os.environ)
            d['L1DSIZE'] = str(val)
            result = subprocess.Popen(run_cmd, env=d, stdout=subprocess.PIPE)
            
            result_out, err = result.communicate()
            result_out = result_out.decode('utf-8')
            
            print(result_out)

            ms_time = float(ms_time_re.search(result_out)[0].split('ms')[0])

            if(times.get(program)==None):
                    times[program] = {val:[[size], [ms_time]]}
            else:
                if(times.get(program).get(val) != None):
                    times[program][val][0].append(size)
                    times[program][val][1].append(ms_time)
                else:
                    times[program][val] = [[size], [ms_time]]

        lin_dic = {'size':times[program][val][0], 'time':times[program][val][1]}
        ax = sns.lineplot(x='size', y='time', legend='full', label=val, markers=True, marker="o", data=lin_dic)
        ax.set(xlabel='size', ylabel='time')


    sns.set(style='darkgrid')
        
    name = "times_for_blocksize.png"
    plt.legend()
    plt.savefig(name)
    plt.clf()

def valgrind_calls_cpp():

    times = {}
    program = 'cpp'
    vals = [128, 256, 512]

    L1, L2 = cache_line()

    double_size = 8

    vals += [int(L1*1024//double_size)]
    vals += [int((L1*1024//double_size)**0.5)]

    for val in vals:
        for size in [128, 256]:
            name_file = 'valgrind_block_s'+str(size)+'_val'+str(val)


            env_var = 'L1DSIZE='+str(val)
            run_cmd = ['valgrind', '--log-file='+name_file, '--tool=cachegrind', './matrix.out', 'block', str(size)] #valgrind --tool=cachegrind

            d = dict(os.environ)
            d['L1DSIZE'] = str(val)
            result = subprocess.Popen(run_cmd, env=d, stdout=subprocess.PIPE)
            
            result_out, err = result.communicate()
            result_out = result_out.decode('utf-8')


if __name__ == '__main__':
    print('Running all CPAR tests and getting results.')
    all_times() # get all running times of different implementations and sizes
    plot_graph(times) # save all the graphs
    cache_sizes_block() # block algorithm with cache size differences
    valgrind_calls_cpp() # save the output of the valgrind valgrind_calls_cpp