import sys
from subprocess import Popen, PIPE
from os import listdir
from os.path import isfile, join
import time
import numpy as np

if len(sys.argv) < 2:
    print("USAGE - python3 tester.py [number_of_runs_per_tests] [optional - tests numbers seperated by space]")

number_of_runs = int(sys.argv[1])
cmds = sys.argv[2:]


dir = 'tests'
print("run tests")
tests = [join(dir, f) for f in listdir(dir) if isfile(join(dir, f)) and join(dir, f).startswith(dir + '/test')]

for test in tests:
    if len(cmds) != 0 and not test.split('tests/test')[1] in cmds:
        continue

    print('running %s' % str(test))

    runs_time = []
    for _ in range(number_of_runs):

        myinput = open(test)
        p = Popen(['./RodPathFinder', 'output'], stdin=myinput, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate()
        time.sleep(0.1)

        if p.returncode != 0:
            print('FAILED %s' % str(err.decode()))
            print('error code %d' % p.returncode)
            continue

        lines = output.decode().split('\n')
        for line in lines:
            if 'Path creation time:' in line:
                runs_time.append(float(line.split('Path creation time:')[1]))
                #print(runs_time[-1])
            if 'Path verifying:' in line:
                if 'SUCCESS' in line:
                    break
        else:
            print('FAILED')
            print(output.decode())
            exit()

    print('%d runs, mean time %f, max is %f, variance %f' %
          (number_of_runs, sum(runs_time)/number_of_runs,max(runs_time), np.var(runs_time)))