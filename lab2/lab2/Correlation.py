import math
import FFT
ACCURACY=8

def MatrixCorrelation(values1,values2,num):
    conMassive = [complex(0,0)]*num
    for m in range(0,num,1):
        for i in range(0,num,1):
            conMassive[m]+=values2[i]*values1[int((m+i) % num)]/num
        conMassive[m] = round(conMassive[m].real,ACCURACY)+1j*round(conMassive[m].imag,ACCURACY)
    return conMassive

def CorrelationWithFT(values1,values2,num):
    w = FFT.GetWForFFT(num)
    massiveY = FFT.FastFT(values1,w,num)
    massiveZ = FFT.FastFT(values2,w,num)
    for i in range(0,num,1):
        massiveZ[i]=massiveZ[i].real + 1j*massiveZ[i].imag*(-1)
    conMassive = []
    for i in range(0,num,1):
        conMassive.append(massiveY[i]*massiveZ[i])
    w = FFT.GetWForOFFT(num)
    resultMassive = FFT.FastFT(conMassive,w,num)
    for i in range(0,num,1):
        resultMassive[i]=round(resultMassive[i].real/(num*num),ACCURACY)+1j*round(resultMassive[i].imag/(num*num),ACCURACY)
    return resultMassive
