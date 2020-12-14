import os
from circuitSimParts import *
from circuitSim import funSPICE

solverOptions = dict(spiceRefNode='0', capAdmittance=0, solverTolerance=10e-6, refV=0.0, maxIters=int(3*10e4))
outputOptions = dict(printRead=False, printResults=False, printSupernodes=False)
codeResults = dict()
expectedResults = dict()
for filename in os.listdir('./testcircuits'):
    if filename[0] != '.' and filename[-4:] != '.asc':
        if 'out' in filename:
            name = filename[:-8]
            expectedResults[name] = readOutput('./testcircuits/' + filename)
        else:
            name = filename[:-4]
            n, e, _ = funSPICE('./testcircuits/' + filename, solverOptions, outputOptions)
            writeOutput(n, e, './output.txt')
            codeResults[name] = readOutput('./output.txt')

for key in codeResults:
    print(key)
    print('     nodes:')
    for k, v in compareOutputs(codeResults[key], expectedResults[key])[0].items():
        if v > 0.01:
            print(k + ': ' + str(v))
    print('     elems:')
    for k, v in compareOutputs(codeResults[key], expectedResults[key])[1].items():
        if v > 0.01:
            print(k + ': ' + str(v))
