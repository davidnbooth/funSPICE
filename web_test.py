import circuitSim
import circuitSimParts
import circuitSimHelpers


input = "R2 N003 0 3\nD1 N002 N003 D\nV1 N001 N004 7\nR1 N002 N001 2\nR3 0 N004 1"
input = 'R2 N003 0 3\r\nD1 N002 N003 D\r\nV1 N001 N004 7\r\nR1 N002 N001 2\r\nR3 0 N004 1'

solverOptions = dict(spiceRefNode='0',
                     capAdmittance=0,
                     solverTolerance=10e-6,
                     refV=0.0,
                     maxIters=int(3*10e4),
                     wrelax=1)
outputOptions = dict(printRead=True,
                     printResults=True,
                     printSupernodes=True)

nodeDict, elemDict, shortedElems = circuitSim.funSPICE(input, solverOptions, outputOptions)
output_text = circuitSimParts.writeOutput(nodeDict, elemDict, '', file_write=False)
print(output_text)
