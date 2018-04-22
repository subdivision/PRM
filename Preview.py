import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import Slider
from matplotlib import animation

obs_lines = []


def draw_obstacles(obstacles, ax):
    for pol in obstacles:
        line1, = ax.plot([], [], 'brown', lw=2)
        obs_lines.append(line1)


if len(sys.argv) < 2:
    print("[USAGE1] output.txt")
    exit()

cmds = sys.argv[1:]
with open(cmds[0], 'r') as myfile:
    data = myfile.read().replace('\n', ' ')
    data = data.replace('\r', ' ')
    data = data.replace('\t', ' ')
    data = data.replace('  ', ' ')
    data = data.replace('  ', ' ')
s = data.split()
i = 0
obs_num = int(s[i])
i += 1
obs = []
for _ in range(obs_num):
    ver_num = int(s[i])
    i += 1
    poly = []
    for __ in range(ver_num):
        poly.append((float(s[i]), float(s[i + 1])))
        i += 2
    poly.append(poly[0])
    obs.append(poly)
rodLength = float(s[i])
i += 1
path = []
pathLength = int(s[i])
i += 1
for _ in range(pathLength):
    path.append(((float(s[i]), float(s[i + 1])), float(s[i + 2]), int(s[i + 3])))
    i += 4

fig = plt.figure()

ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-100, 100), ylim=(-100, 100))

draw_obstacles(obs, ax)
path_num_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
line, = ax.plot([], [], 'o-', lw=2)
frames_per_path_part = 3000


def init():
    """initialize animation"""
    global obs, obs_lines
    ret = [line, path_num_text]
    line.set_data([], [])
    i = 0
    for pol in obs:
        xs, ys = zip(*pol)  # create lists of x and y values
        obs_lines[i].set_data(xs, ys)
        ret.append(obs_lines[i])
        i += 1

    path_num_text.set_text('')
    return tuple(ret)


currentFrame = 0
currentSpeed = 20


def animate(i):
    """perform animation step"""
    global path, frames_per_path_part, rodLength, obs_lines, currentFrame, currentSpeed
    currentFrame = (currentFrame + currentSpeed) % (frames_per_path_part * (pathLength - 1))
    index = int(currentFrame / frames_per_path_part)
    part_of = currentFrame % frames_per_path_part
    part_of = part_of / frames_per_path_part
    p1, r1, _ = path[index]
    p2, r2, cc = path[index + 1]
    r1 = math.fmod(r1, 2 * math.pi)
    r2 = math.fmod(r2, 2 * math.pi)
    if r1 < 0:
        r1 += 2 * math.pi
    if r2 < 0:
        r2 += 2 * math.pi
    if cc == 1:
        if r2 < r1:
            r2 += 2 * math.pi
    else:
        if r2 > r1:
            r1 += 2 * math.pi
    r3 = r1 + (r2 - r1) * part_of
    r3 = math.fmod(r3, 2 * math.pi)

    x = np.cumsum([p1[0] + (p2[0] - p1[0]) * part_of,
                   rodLength * math.cos(r3)], dtype='float')
    y = np.cumsum([p1[1] + (p2[1] - p1[1]) * part_of,
                   rodLength * math.sin(r3)], dtype='float')

    line.set_data(x, y)
    ret = [line, path_num_text]

    for pol in obs_lines:
        ret.append(pol)

    path_num_text.set_text('path part = %d -> %d' % (index, index + 1))
    return tuple(ret)


ani = animation.FuncAnimation(fig, animate, frames=frames_per_path_part * (pathLength - 1), interval=33, blit=True,
                              init_func=init)

axcolor = 'lightgoldenrodyellow'
plt.subplots_adjust(bottom=0.2)
axPlayer = plt.axes([0.1, 0.1, 0.8, 0.025], facecolor=axcolor)
axSpeed = plt.axes([0.1, 0.05, 0.8, 0.025], facecolor=axcolor)


def sliderChange(val):
    global currentFrame
    currentFrame = val


def sliderChangeSpeed(val):
    global currentSpeed
    currentSpeed = val


splayer = Slider(axPlayer, 'Player', 0, frames_per_path_part * (pathLength - 1), valinit=currentFrame)
splayerSpeed = Slider(axSpeed, 'Speed', 0, 100, valinit=currentSpeed)
splayer.on_changed(sliderChange)
splayerSpeed.on_changed(sliderChangeSpeed)

plt.show()
