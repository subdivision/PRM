import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

obs_lines = []

class cPoint():
    def __init__(self, str, rodLength):
        data = str.split()
        self.startX = float(data[0])
        self.startY = float(data[1])
        dir = float(data[2])
        self.endX = self.startX + math.cos(dir) * rodLength
        self.endY = self.startY + math.sin(dir) * rodLength

    def draw(self, color='green'):
        plt.plot([self.startX, self.endX], [self.startY, self.endY], color)


def draw_obstacles(obstacles):
    for pol in obstacles:
        xs, ys = zip(*pol)  # create lists of x and y values
        plt.plot(xs, ys, 'brown')

def load_obs(data):
    res = []
    for obs_str in data[1:]:
        if len(obs_str.strip()) == 0:
            continue
        obs_data = obs_str.split()
        pol = []
        for x,y in zip(obs_data[1::2], obs_data[2::2]):
            pol.append([float(x), float(y)])

        pol.append(pol[0])  # repeat the first point to create a 'closed loop'
        res.append(pol)
    return res

if len(sys.argv) < 2:
    print("[USAGE1] output.txt")
    exit()


cmds = sys.argv[1:]
present_level = 1
if len(cmds) > 1:
    present_level = int(cmds[1])

plt.figure()


data = open(cmds[0], 'r').read()

lines = data.split('\n')
rodLength = float(lines[0])
startCpoint = cPoint(lines[1], rodLength)
endCpoint = cPoint(lines[2], rodLength)

obs = load_obs(lines[4:])


data = open('points', 'r').read()
lines = data.split('\n')
colors = ['black', 'gray', 'orange', 'green']
lines_types = [[],[],[],[]]
for line in lines:
    if len(line.strip()) == 0:
        continue
    type = int(line[0])
    cp = cPoint(line[2:], rodLength)
    lines_types[type].append(cp)

for i, line_type in enumerate(lines_types):
    if i<present_level:
        continue

    for cp in line_type:
        cp.draw(colors[i])

draw_obstacles(obs)
startCpoint.draw(color='blue')
endCpoint.draw(color='red')

plt.show()
