import os


folder_root = "/home/glenn/cluster/benchmark5_attract/*"
file_copy = "*dock.result"

benchmarks =[]
benchmarks.append(('benchmark_GPU_scorig_50cut_0modes','benchmark_GPU_scGPU_0modes'))
benchmarks.append(('benchmark_GPU_scorig_50cut_5modes','benchmark_GPU_scGPU_5modes'))
benchmarks.append(('benchmark_GPU_scorig_50cut_3modes','benchmark_GPU_scGPU_3modes'))
benchmarks.append(('benchmark_GPU_scorig_50cut_5modes','benchmark_GPU_scGPU_5modes_scnoEV'))
benchmarks.append(('benchmark_GPU_scorig_50cut_3modes','benchmark_GPU_scGPU_3modes_scnoEV'))
benchmarks.append(('benchmark_GPU_scorig_50cut_5modes','benchmark_GPU_scGPU_5modes_sc2EV'))
benchmarks.append(('benchmark_GPU_scorig_50cut_3modes','benchmark_GPU_scGPU_3modes_sc2EV'))
#benchmarks.append(('benchmark_GPU_scorig_50cut_5modes_2EV','benchmark_GPU_5modes_2EV_scGPU_50c_0EV'))

for name_oldBenchmark ,name_newBenchmark in benchmarks:
    bash_command = "for dir in $(ls -d %s ); do   mkdir ${dir}/%s; done"%(folder_root,name_newBenchmark)
    #bash_command = "for dir in $(ls -d %s ); do   echo ${dir}/%s; done"%(folder_root,name_newBenchmark)
    os.system(bash_command)
    bash_command = "for dir in $(ls -d %s ); do   cp ${dir}/%s/%s  ${dir}/%s; done"%(folder_root,name_oldBenchmark,file_copy,name_newBenchmark)
    os.system(bash_command)