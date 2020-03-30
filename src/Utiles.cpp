#include "Utiles.h"
#include "Corps.h"
#include "Systeme.h"
#include "Constantes.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>

using namespace std;



void afficher(const vector<double>& v ){for(int i =0 ; i < (int) v.size() ; i++){cout << i << ") " << v[i] << endl;}}

vector<double> operator-(vector<double> a,vector<double>b){
	
	vector<double> difference(3);

	for(int i=0;i<3;i++)
	{
		difference[i]=a[i]-b[i];
	}	

	return difference;
}

double norme(const vector<double>& v){
	double nor(0);
	for(int i = 0 ; i < (int) v.size() ; i++){
		nor += v[i]*v[i];
	}
	nor = sqrt(nor);

	return nor;
}


vector<double> normaliser(const vector<double>& v){
	vector<double> a;	
	a = multScalaire(1.0/norme(v), v);
	return a;
}

double produitScalaire(const vector<double>& a, const vector<double>& b){
	double s(0);	
	for(int i = 0 ; i < (int) a.size() ; i++){
		s+= a[i]*b[i];
	}
	return s;
}

double moyenne(const vector<double>& v){
	double moy;
	for(int i = 0 ; i < (int) v.size() ; i++){
		moy += v[i];
	}
	moy /= v.size();
	return moy;
}

vector<double> multScalaire(double m, const vector<double>& a){
	vector <double> v;
	for(int i = 0 ; i < (int) a.size() ; i++){
		v.push_back(a[i]*m);
	}
	return v;
}

vector<double> distance(const vector<double>& a, const vector<double>& b){
	vector<double> distance;
	if(a.size() == b.size()){
		for(int i = 0 ; i < (int) a.size() ; i++){
				distance.push_back(b[i]-a[i]);
		}
	}
	return distance;
}


vector<double> add(const vector<double>& a, const vector<double>& b){
	vector<double> v;
	if(a.size() == b.size()){
		for(int i = 0 ; i < (int) a.size() ; i++){
				v.push_back(b[i]+a[i]);
		}		
	}
	return v;
}

vector<double> oppose(const vector<double>& a){

	vector<double> oppose;
	oppose = multScalaire(-1, a);
	return oppose;
	
}

vector<double> ProdVec(const vector<double>& U,const vector<double>& V){	
	vector<double> V3(U.size());
	
	V3[0]=U[1]*V[2]-U[2]*V[1];
	V3[1]=U[2]*V[0]-U[0]*V[2];
	V3[2]=U[0]*V[1]-U[1]*V[0];

	return V3;
}