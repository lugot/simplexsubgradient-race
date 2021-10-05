from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

blue = (0, 191, 255)
blue = (0.0, 0.749, 1.0)

green = (87, 205, 128)
green = (0.341, 0.804, 0.5)

red = (205, 92, 92)
red = (0.804, 0.361, 0.361)

orange = (255, 127, 80)
orange = (1.0, 0.5, 0.314)

def delete_ext(filename):
    i = 0
    while filename[i] != '.': 
        i+=1

    return filename[:i]

path = 'results'
instances = [delete_ext(f) for f in listdir(path) if isfile(join(path, f))]
instances = list(set(instances))

small_instances_print = []

for instance in instances:
    subgr = pd.read_csv("results/" + instance + ".subgr_log")

    lines = []
    with open("results/" + instance + ".cplex_log") as f:
        lines = f.readlines()

    cplex_dual = []

    i = 0
    for line in lines:
        if len(line) < 9:
            i += 1
            continue
        if line[:9] == "Iteration":
            break
        i += 1

    i += 1

    extra_lines = 0 
    for line in lines[i:-2]:
        line = line.split()
        if len(line) < 4:
            extra_lines += 1
            continue
        try:
            cplex_dual.append(float(line[1]))
        except ValueError:
            try:
                cplex_dual.append(float(line[2]))
            except IndexError:
                print(line)
            except ValueError:
                extra_lines += 1
        

    assert(len(cplex_dual) ==  len(lines[i:-2]) - extra_lines)
    cplex_timesolve = int(lines[-2].split()[1])

    cplex_noiter = False
    if len(cplex_dual) == 0:
        cplex_noiter = True
        print("cplex NO ITERATION")
        continue

    if len(cplex_dual) == 0:
        continue
    cplex_dual = np.array(cplex_dual)
    print("len cplex dual", len(cplex_dual))
    cplex_time = np.arange(0, float(len(cplex_dual)))
    cplex_time *= float(cplex_timesolve / len(cplex_dual))

    subgr_noiter = False
    if subgr.shape[0] == 0:
        subgr_noiter = True
        print("subgradient NO ITERATION")
        continue

    subgr_dual = subgr["dual"]
    subgr_time = subgr["time"]

    print(instance, file=sys.stderr)
    plt.figure(figsize=(8,6))
    plt.xlabel("Time [ms]")
    plt.ylabel("Objective Value")
    plt.plot(cplex_time, cplex_dual, label="simplex", linewidth=1.2, color=orange)
    plt.plot(subgr_time, subgr_dual, label="subgradient", linewidth=1.2, color=blue)
    # plt.hlines(576.924915956562, xmin = 0, xmax = 500, colors = green, linestyle="dashed", linewidth=1)
    plt.legend()
    plt.show()
    
