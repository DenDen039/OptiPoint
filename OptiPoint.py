import os
import numpy as np
from numpy import pi, sin, cos,tan
import pandas as pd
import plotly.express as px
from tkinter.filedialog import askopenfilename
from tkinter import *

def ellipse(a=1, b =1,x_center=0, y_center=0, ax1 = [1, 0],  ax2 = [0,1],  N=300):
    t = np.linspace(0, 2*pi, N)
    xs = a * cos(t)
    ys = b * sin(t)
    R = np.array([ax1, ax2]).T
    xp, yp = np.dot(R, [xs, ys])
    x = xp + x_center 
    y = yp + y_center
    return x, y
def RotateCoord(x,y,angle,x_center,y_center):
    if(0):
        x1 = x*cos(angle)-y*sin(angle)+x_center
        y1 = +x*sin(angle)+y*cos(angle)+y_center
    else:
        x1 = x*cos(angle)+y*sin(angle)+x_center
        y1 = -x*sin(angle)+y*cos(angle)+y_center
    return x1,y1
def PosToGeo(x, y,zone,L0):#перевод в географические координаты
    e2 = 0.006693421623
    e2f = 0.006738525415
    M_PI = 3.14159265358979323846
    a=6378245
    n = float(zone)
    L0 = 6*n-3
    L0 = L0*M_PI/180;
    y1 = y-n*1000000-500000
    x1 = x
    m0 = a*(1-e2)
    m2 = 3/2*e2*m0
    m4 = 5/4*e2*m2
    m6 = 7/6*e2*m4
    m8 = 9/8*e2*m6
    q0 = m0+(1/2)*m2+3/8*m4+5/16*m6+35/128*m8
    q2 = m2*(1/2)+(1/2)*m4+15/32*m6+7/16*m8
    q4 = 1/8*m4+3/16*m6+7/32*m8
    q6 = 1/32*m6+1/16*m8
    b0 = x1/q0
    b2 = -(q2/(2*q0))
    b4 = (q4/(4*q0))
    b6 = -(q6/(6*q0))
    d2 = b2*(1/2*(b2**2)-b4-1)
    d4 = b2**2-b4
    d6 = 3*b2*b4-3/2*(b2**3)-b6
    Bx = b0+d2*sin(2*b0)+d4*sin(4*b0)+d6*sin(6*b0)
    m2x = (cos(Bx)**2)*(e2/(1-e2))
    tx = tan(Bx)
    Nx = a/(1-e2*sin(Bx)*sin(Bx))**(1/2)
    A2 = -(1/(2*Nx**2))*(1+m2x)*tx
    A4 = -(A2/(12/Nx**2))*(5+3*tx**2+m2x-9*m2x*tx*tx-4*m2x*m2x)
    A6 =A2/(360*Nx**4)*(61+90*tx**2+45*tx**4+46*m2x-252*m2x*tx**2-90*m2x*tx**4)
    A8 = -A2/(20160*Nx**6)*(1385+3633*tx**2+4095*tx**4+1575*tx**6)
    B1 = 1/(Nx*cos(Bx))
    B3 = -B1/(6*Nx**2)*(1+2*tx**2+m2x)
    B5 = B1/(120*Nx**4)*(5+28*tx**2+24*tx**4+6*m2x+8*m2x*tx**2)
    B7 = -B1/(5040*Nx**6)*(61+662*tx**2+1320*tx**4+720*tx**6)
    B = Bx+A2*y1**2+A4*y1**4+A6*y1**6+A8*y1**8
    l = B1*y1+B3*y1**3+B5*y1**5+B7*y1**7
    L = l+L0
    bf = x/6367558.497
    b = bf + 0.002518465*sin(2*bf)+0.000003700*sin(4*bf)+0.000000007*sin(6*bf)
    Nf=a/(1-e2*sin(b)**2)**(1/2)
    muf = e2f*cos(b)**2 
    b_B = ((y1**2*(1+muf)*tan(b))/(2*Nf**2))*(1-(y1**2/(12*Nf**2*cos(b)**2))*(3-9*muf+2*cos(b)**2*(1+5*muf-2*muf**2))+(y1**4)/(360*Nf**4*cos(b)**4)*(45-18*muf*(5+4*cos(b)**2)+16*cos(b)**4*(1+13*muf)))
    B= b-(b_B)
    B = (B*180)/M_PI
    L = (L*180)/M_PI
    return B,L

################################################################## Считывание данных
dirct = os.path.abspath(os.curdir)
root = Tk()
filename = askopenfilename()
data = pd.read_csv(filename)#Ввод сsv данных
root.destroy()
file = open("Statistic.txt", "w")

for index, row in data.iterrows():
    file.write(str(int(row['Density']))+' ')
    file.write(str(row['Lat'])+' ')
    file.write(str(row['Lon'])+' ')
    file.write('\n')
file.close()

os.system(dirct+'/CppCore.exe')#Запуск мат ядра
with open('outdata.txt') as file:
        inputf = file.read().split(' ')
file = open('FinalOutput.txt', "w")
Mx,My,r,L0,angle,ae,be,zone = inputf[0],inputf[1],inputf[2],inputf[3],inputf[4],inputf[5],inputf[6],inputf[7]
ax = [float(inputf[8]),float(inputf[9])]
bx = [float(inputf[10]),float(inputf[11])]
if int(inputf[12]) == 0:
    ax,bx=bx,ax
os.remove(dirct+'/outdata.txt')#удаление временных файлов
os.remove(dirct+'/Statistic.txt')#удаление временных файлов
####################################################################################### 

fig = px.scatter_mapbox(data,lat="Lat", lon="Lon", zoom=13,color = "Density")#Создание карты с данными о жителях

################################################################## Создание эллипса
x,y = ellipse(float(ae), float(be))
for i in range(len(x)):
    x[i],y[i] = RotateCoord(x[i],y[i],float(angle),float(Mx),float(My))
    x[i],y[i] = PosToGeo(float(x[i]),float(y[i]),zone,L0)# Перевод координат
fig2 = px.line_mapbox(lat=x, lon=y)
################################################################## 

################################################################# Добавление оптимальной точки
Mx,My = PosToGeo(float(Mx),float(My),zone,L0)
markerx = [Mx]
markery = [My]
markername = ["Optimal point"]
fig3 = px.scatter_mapbox(lon=markery, lat=markerx, size = markery,hover_name = markername)
##########################################################

################################################################## Вывод параметров в тхт файл
file.write('Optimal point: Lat: '+str(Mx)+'; Lon:'+str(My)+'\n')
file.write('Ellipse parameters: '+ ' a = '+str(ae)+ '; b = '+str(be)+'; vector1 = '+str(ax)+'; vector2 = '+str(bx))
file.close()
################################################################## 

##############################################Вывод карты
fig.add_trace(fig2.data[0])
fig.add_trace(fig3.data[0])
fig.update_layout(mapbox = {'style': "open-street-map", 'center': {'lat': Mx,  'lon': My}},margin ={'l':0,'t':0,'b':0,'r':0})
fig.show()