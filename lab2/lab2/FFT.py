import math
ACCURACY=8

def BitReverseSwap(values,num):
    flags=[bool(False)]*num
    kol=int(math.log(num,2))
    for valueI in range(0,num,1):
        if flags[valueI]==False:
            index=valueI
            for i in range(0,int(kol / 2),1):
                lBit = int((valueI / int(math.pow(2,kol-i-1))) % 2)
                rBit = int((valueI / int(math.pow(2,i))) % 2)
                index = index + (rBit - lBit) * int(math.pow(2,kol-i-1)) + (lBit - rBit) * int(math.pow(2,i))
            element = values[index]
            values[index] = values[valueI]
            values[valueI] = element
            flags[index] = True
            flags[valueI] = True
    return values

def GetWForFFT(num):
    row=int(math.log(num,2))
    massive=[]
    for i in range(0,row,1):
        massive.append([complex(0,0)]*int(math.pow(2,i+1)))
    for i in range(0,row,1):
        for j in range(0,int(math.pow(2,i+1)),1):
            massive[i][j]=complex(round(math.cos(2*j*math.pi/int(math.pow(2,i+1))),ACCURACY),round(-1*(math.sin(2*j*math.pi/int(math.pow(2,i+1)))),ACCURACY))
    return massive

def GetWForOFFT(num):
    row=int(math.log(num,2))
    massive=[]
    for i in range(0,row,1):
        massive.append([complex(0,0)]*int(math.pow(2,i+1)))
    for i in range(0,row,1):
        for j in range(0,int(math.pow(2,i+1)),1):
            massive[i][j]=complex(round(math.cos(2*j*math.pi/int(math.pow(2,i+1))),ACCURACY),round(math.sin(2*j*math.pi/int(math.pow(2,i+1))),ACCURACY))
    return massive

# Быстрое преобразование Фурье
def FastFT(values,w,num):
    #переупорядовачивание элементов
    values = BitReverseSwap(values,num)
    mas=[complex(0,0)]*num
    for i in range(0,int(math.log(num,2)),1):
        for j in range(0,num,1):
            if int(j % int(pow(2,i+1))) < int(int(pow(2,i+1))/2):
                mas[j]= values[int(pow(2,i+1)/2)+j] * w[i][int(j % int(pow(2,i+1)))] + values[j]
            else:
                mas[j]= values[j] * w[i][int(j % int(pow(2,i+1)))] + values[int(j-1*int(math.pow(2,i+1)/2))]
        for i in range(0,num,1):
            values[i]= round(mas[i].real,ACCURACY) + 1j*round(mas[i].imag,ACCURACY)
    return values           

def DFT(values,num):
    massive = []
    print("\nDFT")
    for m in range(0,num,1):
        elem = complex(0,0)
        for n in range(0,num,1):
            elem+=values[n]* (math.cos(2*math.pi*m*n/num) - 1j* math.sin(2*math.pi*m*n/num))
        massive.append(complex(round(elem.real,ACCURACY),round(elem.imag,ACCURACY)))
    return massive

def ODFT(values,num):
    massive = []
    print("\nODFT")
    for m in range(0,num,1):
        elem = complex(0,0)
        for n in range(0,num,1):
            elem+=values[n]* (math.cos(2*math.pi*m*n/num) + 1j* math.sin(2*math.pi*m*n/num))
        massive.append(complex(round(elem.real/num,ACCURACY),round(elem.imag/num,ACCURACY)))
    return massive
