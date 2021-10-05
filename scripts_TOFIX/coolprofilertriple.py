from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

blue = (0, 191, 255)
blue = (0.0, 0.749, 1.0)

green = (87, 205, 128)
green = (0.341, 0.804, 0.5)

red = (205, 92, 92)
red = (0.804, 0.361, 0.361)

orange = (255, 127, 80)
orange = (1.0, 0.5, 0.314)

df1 = pd.read_csv("new1")
df2 = pd.read_csv("new2")
df3 = pd.read_csv("new3")

subgr_dual1 = df1["dual"]
subgr_dual2 = df2["dual"]
subgr_dual3 = df3["dual"]

plt.figure(figsize=(6,4))
plt.xlabel("Iteration")
plt.ylabel("Objective Value")

# plt.plot(np.arange(0, 299, 1), subgr_dual1, color = green, linewidth = 1.2, label="correction first 100 iters")
plt.plot(subgr_dual2, color = red, linewidth = 1.2, label="no correction")
plt.plot(subgr_dual3, color = blue, linewidth = 1.2, label="correction")

plt.yscale("symlog")

plt.hlines(576.924915956562, xmin = 0, xmax = 300, colors = green, linestyle="dashed", linewidth=1)
plt.legend()
plt.savefig("correction_log.pdf")

# def delete_ext(filename):
#     i = 0
#     while filename[i] != '.': 
#         i+=1
#
#     return filename[:i]
#
# path = 'results'
# instances = [delete_ext(f) for f in listdir(path) if isfile(join(path, f))]
# instances = list(set(instances))
#
# for instance in instances[:1]:
#     subgr = pd.read_csv("results/" + instance + ".subgr_log")
#     print(subgr)
#
#     lines = []
#     with open("results/" + instance + ".cplex_log") as f:
#         lines = f.readlines()
#
#     cplex_dual = []
#
#     i = 0
#     for line in lines:
#         if len(line) < 9:
#             i += 1
#             continue
#         if line[:9] == "Iteration":
#             break
#         i += 1
#
#     i += 1
#
#     extra_lines = 0 
#     for line in lines[i:-2]:
#         line = line.split()
#         if len(line) < 4:
#             extra_lines += 1
#             continue
#         try:
#             cplex_dual.append(float(line[1]))
#         except ValueError:
#             try:
#                 cplex_dual.append(float(line[2]))
#             except IndexError:
#                 print(line)
#             except ValueError:
#                 extra_lines += 1
#         
#
#     assert(len(cplex_dual) ==  len(lines[i:-2]) - extra_lines)
#     cplex_timesolve = int(lines[-2].split()[1])
#
#     if len(cplex_dual) == 0:
#         print("cplex NO ITERATION")
#         continue
#
#     cplex_dual = np.array(cplex_dual)
#     cplex_time = np.linspace(start=0, stop=cplex_timesolve, num=len(cplex_dual))
#
#     if subgr.shape[0] == 0:
#         print("subgradient NO ITERATION")
#         continue
#
