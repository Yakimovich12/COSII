import math
import FFT
import pylab
from matplotlib import mlab
import matplotlib
import Convolution as Con
import Correlation as Cor
ACCURACY=8

def CopyMassive(massive,num):
    resultMas=[]
    for i in range(0,num,1):
        resultMas.append(massive[i])
    return resultMas


def GetFunction1Values(num):
    values = list()
    for i in range(0,num,1):
        values.append(complex(round(math.sin(2*math.pi*i/360),ACCURACY)))
    return values

def GetFunction2Values(num):
    values = list()
    for i in range(0,num,1):
        values.append(complex(round(math.cos(2*math.pi*4*i/360),ACCURACY)))
    return values

def DrawGraphicFunction():
    xVal=[]*360
    yVal1=[]*360
    yVal2=[]*360
    for i in range(0,360,1):
        xVal.append(2*math.pi*i/360)
        yVal1.append(math.sin(xVal[i]))
        yVal2.append(math.cos(4*xVal[i]))
    pylab.figure(0)
    pylab.plot(xVal,yVal1,label="f(x)=sin(x)")
    pylab.plot(xVal,yVal2,label="f(x)=cos(4x)")
    pylab.text(2,1.2,u"Графики исходных функций")
    legend=pylab.legend()
    legend.set_draggable(True)
    pylab.grid()
    pylab.show()

def DrawValues(values,num,message,figure):
    yVal=[]
    xVal=[]
    for i in range(0,num,1):
        yVal.append(math.sqrt(math.pow(values[i].real,2)+math.pow(values[i].imag,2)))
        xVal.append(i*2*math.pi/360)
    pylab.figure(figure)
    pylab.plot(xVal,yVal,label=message)
    legend=pylab.legend()
    legend.set_draggable(True)
    pylab.grid()
    pylab.show()



def FT(values1,values2,num):
     w=FFT.GetWForFFT(num)
     massive1=FFT.FastFT(values1,w,num)
     massive2=FFT.FastFT(values2,w,num)

     DrawValues(massive1,num,"Дискретные значения после БПФ для f(x)=sin(x)",1)
     print("\nДискретные значения после быстрого преобразования Фурье для f(x)=sin(x)")
     for i in range(0,num,1):
         print("Действительная часть: ",massive1[i].real," Мнимая часть: ",massive1[i].imag)

     DrawValues(massive2,num,"Дискретные значения после БПФ для f(x)=cos(4x)",1)
     print("\nДискретные значения после быстрого преобразования Фурье для f(x)=cos(4x)")
     for i in range(0,num,1):
         print("Действительная часть: ",massive2[i].real," Мнимая часть: ",massive2[i].imag)

     w=FFT.GetWForOFFT(num)
     result1=FFT.FastFT(massive1,w,num)
     result2=FFT.FastFT(massive2,w,num)
     for i in range(0,num,1):
         result1[i]=round(result1[i].real/num,ACCURACY)+1j*round(result1[i].imag/num,ACCURACY)
     for i in range(0,num,1):
         result2[i]=round(result2[i].real/num,ACCURACY)+1j*round(result2[i].imag/num,ACCURACY)

     DrawValues(result1,num,"Дискретные значения после ОБПФ (исходные значения для f(x)=sin(x))",2)
     print("\nДискретные значения сигнала после обратного быстрого преобразования Фурье (исходные данные для f(x)=sin(x))")
     for i in range(0,num,1):
         print("Действительная часть: ",result1[i].real," Мнимая часть: ",result1[i].imag)

     DrawValues(result2,num,"Дискретные значения после ОБПФ (исходные значения для f(x)=cos(4x))",2)
     print("\nДискретные значения сигнала после обратного быстрого преобразования Фурье (исходные данные для f(x)=cos(4x))")
     for i in range(0,num,1):
         print("Действительная часть: ",result2[i].real," Мнимая часть: ",result2[i].imag)

def ConvolutionAndCorrelation(values1,values2,num):
    conMat = Con.MatrixConvolution(values1,values2,num)
    mas1=CopyMassive(values1,num)
    mas2=CopyMassive(values2,num)
    conFFT = Con.ConvolutionWithFT(mas1,mas2,num)
    DrawValues(conFFT,num,"Вычисление свертки функций с помощью быстрого преобразования Фурье",3)
    DrawValues(conMat,num,"Вычисление свертки функций с помощью матриц",3)
    print("\nВычисление свертки функций с помощью матриц")
    for i in range(0,num,1):
        print("Действительная часть: ",conMat[i].real," Мнимая часть: ",conMat[i].imag)
    print("\nВычисление свертки функций с помощью быстрого преобразования Фурье")
    for i in range(0,num,1):
        print("Действительная часть: ",conFFT[i].real," Мнимая часть: ",conFFT[i].imag)

    corMat = Cor.MatrixCorrelation(values1,values2,num)
    mas1=CopyMassive(values1,num)
    mas2=CopyMassive(values2,num)
    corFFT = Cor.CorrelationWithFT(mas1,mas2,num)
    DrawValues(corFFT,num,"Вычисление корреляции функций с помощью быстрого преобразования Фурье",4)
    DrawValues(corMat,num,"Вычисление корреляции функций с помощью матриц",4)
    print("\nВычисление корреляции функций с помощью матриц")
    for i in range(0,num,1):
        print("Действительная часть: ",corMat[i].real," Мнимая часть: ",corMat[i].imag)
    print("\nВычисление корреляции функций с помощью быстрого преобразования Фурье")
    for i in range(0,num,1):
        print("Действительная часть: ",corFFT[i].real," Мнимая часть: ",corFFT[i].imag)

def main():
    print("Введите количество отсчетов")
    num=int(input())
    # Получаем отсчеты для данной функции
    values = GetFunction1Values(num)
    values2 = GetFunction2Values(num)
    # Выводим отсчеты заданной функции
    print("\nВыводим отсчеты заданной функции 1")
    for i in range(0,num,1):
        print(values[i].real)
    print("\nВыводим отсчеты заданной функции 2")
    for i in range(0,num,1):
        print(values2[i].real)
    DrawGraphicFunction()
    val1=CopyMassive(values,num)
    val2=CopyMassive(values2,num)
    FT(val1,val2,num)
    ConvolutionAndCorrelation(values,values2,num)
    

main()

