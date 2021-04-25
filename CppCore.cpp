#include <bits/stdc++.h>
using namespace std;
int Data_Size = 1129;//количество точек
#define M_PI 3.14159265358979323846//Pi
struct Point{//структура точка
	double x,y;
};
struct PointPeople{//структура точка с частотой
	Point coord;
	double freq;
	void PrintData(){
		cout << "|" << freq;
		cout << "|";
		printf("%.12f",coord.x);
		cout << "|";
		printf("%.12f",coord.y);
		cout << "|" << endl;
	}
};
PointPeople data[100000] = {};//массив всех точек с частотой
double l0 = 42,zone = 7;
double* equation(double a, double b, double c, double d) {//квадратное уравнение вида: (a-x)*(b-x)-c*d
    double z = -1 * (b + a);
    double v = a * b - c * d;
    double D = z * z - 4 * v;
    double* x;
    x = new double[2];
    if (D < 0) {
        cout << "Error";
        return 0;
    }
    else {
        x[0] = (-1 * z + sqrt(z * z - 4 * v)) / 2; 
        x[1] = (-1 * z - sqrt(z * z - 4 * v)) / 2; 
        return x;
    }
}

double* SLAU(double** a, double x2) {//Решение слау параметры массив векторов и свободная переменная
    double* x;
    x = new double[2];
    x = equation(a[0][0], a[1][1], a[1][0], a[0][1]);
    double* ans;
    ans = new double[4];
    int j = 0;
    for (int i = 0; i < 4; i += 2) {        
        double x1 = -a[0][1] / (a[0][0] - x[j]);
        ans[i] = x1;
        ans[i + 1] = x2;
        j++;
    }
    return ans;
}

double dot_product(double *x, double *y) {//скалярное произведение параметры вектор1 вектор2
    int i; 
    double ans = 0; 
    for(i=0; i<2; ++i) 
        ans += x[i]*y[i]; 
    return ans; 
} 
void normalize(double *x) { //нормирование вектора параметры вектор
    double norm = sqrt(dot_product(x, x)); 
    int i; 
    for(i=0; i<2; ++i) 
        x[i] /= norm; 
} 
void gram_schimdt(double q[][2], int n) { 
    int i, j, k; 
    for(i=1; i<n; ++i) { 
        for(j=0; j<i; ++j) { 
            double scaling_factor = dot_product(q[j], q[i])/ dot_product(q[j], q[j]); 
            for(k=0; k<2; ++k) 
                q[i][k] -= scaling_factor*q[j][k]; 
        } 
    } 
    for(i=0; i<n; ++i) 
        normalize(q[i]); 
    } 
     

double* FindVectors(int a, int b, int c, int d) {//Поиск напр. векторов матрицы
    double** arr = new double* [2];
    for (int i = 0; i < 2; i++) {
        arr[i] = new double[2];
    }
    arr[0][0] =a;
    arr[0][1] =b;
    arr[1][0] =c;
    arr[1][1] =d;
    double* s = SLAU(arr, 1);
   	delete arr;
    return s;
    
}
Point GeoToPos(double B, double L){//переобразование гео координат в прямоугольные. Параметры широта и долгота
	long double  n = zone;
	B = B * M_PI / 180;
	L = L * M_PI / 180;
  	long double e2 = 0.006693421623;//2 ексцентриситет
  	long double e2f = 0.006738525415;//1 ексцентриситет
  	long double a = 6378245;//большая ось эллипсоида
  	long double mu = e2f * cos(2 * B);//мю
  	long double N = a / (sqrt(1 - e2 * sin(B) * sin(B)));
  	long double X = 6367558.497 * (B - 0.002518466 * sin(2 * B) + 0.000002643 * sin(4 * B) - 0.000000003 * sin(6 * B));
  	long double L0 = 6 * n - 3;
  	L0 = L0 * M_PI / 180;
  	long double l = L - L0;
  	long double x = X + 1 / 4.0 * N * l * l * sin(2 * B) * (1 - 1 / 12.0 * l * l * (1 - cos(B) * cos(B) * (6 + 9 * mu + 4 * mu * mu)) + 1 / 360.0 * pow(l, 4) * (1 - 60 * cos(B) * cos(B) + 120 * pow(cos(B), 4) - 330 * mu * cos(B) * cos(B) + 600 * mu * cos(B)));
  	long double y = N * l * cos(B) * (1 - 1 / 6.0 * l * l * (1 - (2 + mu) * cos(B) * cos(B)) + 1 / 120.0 * pow(l, 4) * (1 - (20 - 24 * cos(B) * cos(B) + 58 * mu - 72 * mu * cos(B) * cos(B)) * cos(B) * cos(B)));
  	y = n * 1000000 + 500000 + y;
	Point p;
	p.x = x; p.y = y;
	return p;
}
void InputData(){//ввод данных
	ifstream stat("statistic.txt");//база данных=)
	double lat,lon;//долгота и широта
	Data_Size = 0;
	while(true){
		if(stat.eof())
   			break;
		stat >>	data[Data_Size].freq;
		stat >> lat;
		stat >> lon;
		if(l0 == -1){
			l0 = floor(lon)+1;
			while(int(l0) % 6 != 0)
				l0++;
			zone=int(l0)/6;
		}
			
		data[Data_Size].coord = GeoToPos(lat,lon);
		Data_Size++;
	}
	
}
Point Dispersion (PointPeople *arr, int N,Point expect){//дисперсия xy параметры массив точек с частотой и размер массива
	Point desp;
	desp.x = 0;
	desp.y = 0;
	double people = 0;//общие значение людей
	for(int i = 0; i < N;i++){
		people+=arr[i].freq;
	}
	for(int i = 0; i < N;i++){
		desp.x = desp.x+(arr[i].coord.x-expect.x)*(arr[i].coord.x-expect.x)*(arr[i].freq/people);
		desp.y = desp.y+(arr[i].coord.y-expect.y)*(arr[i].coord.y-expect.y)*(arr[i].freq/people);
	}
	return desp;
}
double Expectationxy(PointPeople *arr, int N){//мат ожидание xy параметры массив точек с частотой и размер массива
	double people = 0;//общие значение людей
	for(int i = 0; i < N;i++){
		people+=arr[i].freq;
	}
	double expect = 0;
	for(int i = 0; i < N;i++){
		expect = expect + arr[i].coord.y*arr[i].coord.x*(arr[i].freq/people);
	}
	return expect;
}
Point Expectations(PointPeople *arr, int N){//мат ожидание параметры массив точек с частотой и размер массива
	double people = 0;//общие значение людей
	for(int i = 0; i < N;i++){
		people+=arr[i].freq;
	}
	Point expect;
	expect.x = 0;
	expect.y = 0;
	for(int i = 0; i < N;i++){
		expect.x = expect.x+ arr[i].coord.x*arr[i].freq/people;
		expect.y = expect.y+ arr[i].coord.y*arr[i].freq/people;
	}
	return expect;
}
Point RotateCoord(Point p, double alpha){//поворот координат параметры точка и угол
	Point pt;
	pt.x = p.x*cos(alpha)+p.y*cos(alpha);
	pt.y =-p.x*cos(alpha)+p.y*cos(alpha);
	return pt;
}
double Angle (double *v1, double *v2){//угол между векторами 
	return dot_product(v1,v2)/(sqrt(dot_product(v1,v1))*sqrt(dot_product(v2,v2)));
}
int main (){
	ofstream out("outdata.txt");
	InputData();//ввод данных
	Point M = Expectations(data,Data_Size);//мат ожидание
	Point D = Dispersion (data,Data_Size,M);//дисперсия 
	double Mxy = Expectationxy(data,Data_Size);//мат ожидание xy
	double r = (Mxy-M.x*M.y)/(sqrt(D.x)*sqrt(D.y));//корреляция
	double alpha = atan((2*r*sqrt(D.x)*sqrt(D.y))/(D.x*D.x-D.y*D.y))/2;//угол поворота оси
	double ae = sqrt(2/(1/D.x+1/D.y-sqrt(pow((1/D.x-1/D.y),2)+4*(r/(sqrt(D.x)*sqrt(D.y))))*(r/(sqrt(D.x)*sqrt(D.y)))));//большая ось элипса рассеивания
	double be = sqrt(2/(1/D.x+1/D.y+sqrt(pow((1/D.x-1/D.y),2)+4*(r/(sqrt(D.x)*sqrt(D.y))))*(r/(sqrt(D.x)*sqrt(D.y)))));//малая ось элипса рассеивания
	double *arr = FindVectors(D.x,r*sqrt(D.x)*sqrt(D.y),r*sqrt(D.x)*sqrt(D.y),D.y);
	
	/////////////////////нахождение векторов
	double matrix[][2]={{arr[0],arr[1]},{arr[2],arr[3]}};
	double iv[2] = {1,0};
	double v1[2] ={arr[0],arr[1]};
	double v2[2] ={arr[2],arr[3]};
    gram_schimdt(matrix,2);
    //////////////////////
    
	out << M.x << " " << M.y << " " << r << " " << l0 << " " << alpha << " " << ae << " " << be << " " << zone << " " << matrix[0][0]  << " "<< matrix[0][1]  << " "<< matrix[1][0]   << " "<< matrix[1][1] ;
	if(abs(alpha-Angle(v1,iv)) > abs(alpha-Angle(iv,v2)))
		out << " 1";
	else
		out << " 0";
	return 0;
}

