import sys

def f(filename):
    result = []
    with open(filename) as file:
        for lines in file:
            result += lines.split()
    return map(float, result)



l1 = f(sys.argv[1])
l2 = f(sys.argv[2])

for i in xrange(len(l1)):
    if l1[i]-l2[i] > 0.01:
        exit(1)
