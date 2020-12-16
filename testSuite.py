import os
from circuitSimParts import *
from circuitSim import funSPICE

testTolerance = 0.008
solverOptions = dict(spiceRefNode='0', capAdmittance=0, solverTolerance=10e-6, refV=0.0, maxIters=int(3*10e4), wrelax=1)
outputOptions = dict(printRead=False, printResults=False, printSupernodes=False)
codeResults = dict()
expectedResults = dict()
for filename in os.listdir('./testcircuits'):
    if filename[0] != '.' and filename[-4:] == '.txt':
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
    print('   Node Voltage Errors:')
    print('      %      abs. [V]')
    for k, v in compareOutputs(codeResults[key], expectedResults[key])[0].items():
        if v[0] > testTolerance:
            print(k + ': ' + str(round(v[0]*100, 3)) + ', ' + str(round(v[1], 3)))
    print('   Element Current Errors:')
    print('      %      abs. [V]')
    for k, v in compareOutputs(codeResults[key], expectedResults[key])[1].items():
        if v[0] > testTolerance:
            print(k + ': ' + str(round(v[0]*100, 3)) + ', ' + str(round(v[1], 3)))
