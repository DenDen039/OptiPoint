#include <bits/stdc++.h>
using namespace std;
int Data_Size = 1129;//���������� �����
#define M_PI 3.14159265358979323846//Pi
struct Point{//��������� �����
	double x,y;
};
struct PointPeople{//��������� ����� � ��������
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
PointPeople data[100000] = {};//������ ���� ����� � ��������
double l0 = 42,zone = 7;

double* equation(double a, double b, double c, double d) {//���������� ��������� ����: (a-x)*(b-x)-c*d
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
double* SLAU(double** a, double x2) {//������� ���� ��������� ������ �������� � ��������� ����������
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
double dot_product(double *x, double *y) {//��������� ������������ ��������� ������1 ������2
    int i; 
    double ans = 0; 
    for(i=0; i<2; ++i) 
        ans += x[i]*y[i]; 
    return ans; 
} 
void normalize(double *x) { //������������ ������� ��������� ������
    double norm = sqrt(dot_product(x, x)); 
    int i; 
    for(i=0; i<2; ++i) 
        x[i] /= norm; 
}     
double* FindVectors(double a, double b, double c, double d) {//����� ����. �������� �������
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

Point GeoToPos(double B, double L){//��������������� ��� ��������� � �������������. ��������� ������ � �������
	long double  n = zone;
	B = B * M_PI / 180;
	L = L * M_PI / 180;
  	long double e2 = 0.006693421623;//2 ��������������
  	long double e2f = 0.006738525415;//1 ��������������
  	long double a = 6378245;//������� ��� ����������
  	long double mu = e2f * cos(2 * B);//��
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
void InputData(){//���� ������
	ifstream stat("Statistic.txt");//���� ������=)
	double lat,lon;//������� � ������
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
Point Dispersion (PointPeople *arr, int N,Point expect){//��������� xy ��������� ������ ����� � �������� � ������ �������
	Point desp;
	desp.x = 0;
	desp.y = 0;
	double people = 0;//����� �������� �����
	for(int i = 0; i < N;i++){
		people+=arr[i].freq;
	}
	for(int i = 0; i < N;i++){
		desp.x = desp.x+(arr[i].coord.x-expect.x)*(arr[i].coord.x-expect.x)*(arr[i].freq/people);
		desp.y = desp.y+(arr[i].coord.y-expect.y)*(arr[i].coord.y-expect.y)*(arr[i].freq/people);
	}
	return desp;
}
double Expectationxy(PointPeople *arr, int N){//��� �������� xy ��������� ������ ����� � �������� � ������ �������
	double people = 0;//����� �������� �����
	for(int i = 0; i < N;i++){
		people+=arr[i].freq;
	}
	double expect = 0;
	for(int i = 0; i < N;i++){
		expect = expect + arr[i].coord.y*arr[i].coord.x*(arr[i].freq/people);
	}
	return expect;
}
Point Expectations(PointPeople *arr, int N){//��� �������� ��������� ������ ����� � �������� � ������ �������
	double people = 0;//����� �������� �����
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
double Angle (double *v1, double *v2){//���� ����� ��������� 
	return dot_product(v1,v2)/(sqrt(dot_product(v1,v1))*sqrt(dot_product(v2,v2)));
}
int main (){
	ofstream out("outdata.txt");
	InputData();//���� ������
	Point M = Expectations(data,Data_Size);//��� ��������
	Point D = Dispersion (data,Data_Size,M);//��������� 
	double Mxy = Expectationxy(data,Data_Size);//��� �������� xy
	double r = (Mxy-M.x*M.y)/(sqrt(D.x)*sqrt(D.y));//����������
	double alpha = atan((2*r*sqrt(D.x)*sqrt(D.y))/(D.x-D.y))/2;//���� �������� ���
	long double A = 1/D.x;
	long double B = 1/D.y;
	long double C = r/(sqrt(D.x)*sqrt(D.y));
	double ae = sqrt(2/(A+B-sqrt((A-B)*(A-B)+4*C*C)));//������� ��� ������ �����������
	double be = sqrt(2/(A+B+sqrt((A-B)*(A-B)+4*C*C)));//����� ��� ������ �����������
	/////////////////////���������� ��������
	double *arr = FindVectors(D.x,r*sqrt(D.x)*sqrt(D.y),r*sqrt(D.x)*sqrt(D.y),D.y);
	double iv[2] = {1,0};
	double v1[2] ={arr[0],arr[1]};
	double v2[2] ={arr[2],arr[3]};
	normalize(v1);
	normalize(v2);
    //////////////////////
	out << M.x << " " << M.y << " " << r << " " << l0 << " " << alpha << " " << ae << " " << be << " " << zone << " " << v1[0]  << " "<< v1[1]  << " "<< v2[0]   << " "<< v2[1] ;
	if((abs(abs(alpha)-abs(Angle(iv,v1))) > abs(abs(alpha)-abs(Angle(iv,v2))) && alpha > 0) || (abs(abs(alpha)-abs(Angle(iv,v1))) < abs(abs(alpha)-abs(Angle(iv,v2))) && alpha < 0))
		out << " 1";
	else
		out << " 0";
	return 0;
}

