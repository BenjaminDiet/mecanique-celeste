#include "Corps.h"
#include "Utiles.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;


// Constructeurs
Corps::Corps() : position(3, 0), vitesse(3,0), acc(3,0), masse(0)
{}


Corps::Corps(vector<double> r, vector<double> v, vector<double>a, double m, string s, string n, int nat) : position(r), vitesse(v), acc(a), masse(m), lien(s), nom(n), nature(nat)
{}



// Accesseurs
vector<double> Corps::getPosition() const{return position;}
vector<double> Corps::getVitesse() const{return vitesse;} 
vector<double> Corps::getAcc() const{return acc;}
double Corps::getMasse() const{return masse;}
double Corps::getAire() const{return aire;}
string Corps::getLien() const{return lien;}
int Corps::getNature() const{return nature;}


// Autres méthodes
void Corps::addAcc(vector <double> acce){for(size_t k = 0 ; k < acce.size() ; k++){this->acc[k] += acce[k];}}
void Corps::setAcc(vector <double> acce){for(size_t k = 0 ; k < acce.size() ; k++){this->acc[k] = acce[k];}}
void Corps::saveAcc(){this->accAvant.assign(this->acc.begin(), this->acc.end());}
void Corps::emptyAcc(){for(size_t i = 0 ; i < this->acc.size() ; i++){this->acc[i] = 0;}}


void Corps::loiDesAires(vector <double> barycentre, double h){
	vector<double> momentCinetique(3);
	momentCinetique = ProdVec(distance(this->position, barycentre), this->vitesse);
	this->aire = norme(momentCinetique)*h/2.0;
}


// Résolution


// Verlet
void Corps::majPositionVerlet(double h){
	for(int i = 0 ; i < (int) this->position.size() ; i++){
		this->position[i] += this->vitesse[i]*h +0.5*this->acc[i]*h*h;
	}
}

void Corps::majVitesseVerlet(double h){
	for(int i = 0 ; i < (int) this->vitesse.size() ; i++){
		this->vitesse[i] += 0.5*(this->acc[i]+this->accAvant[i])*h;
	}
}




// Euler
void Corps::majPositionEuler(double h){
	for(int i = 0 ; i < (int) this->position.size() ; i++){
		this->position[i] += this->vitesse[i]*h;
	}
}



void Corps::majVitesseEuler(double h){
	for(int i = 0 ; i < (int) this->vitesse.size() ; i++){
		this->vitesse[i] += this->acc[i]*h; 
	}
}








// RK4
void Corps::AddPosition(vector<double> k,double h,double indice){
		for(int i=0;i<(int)position.size();i++){position[i]=position[i]+k[indice+i]*h;}
}

void Corps::SetPosition(vector<double> k,double indice){
		for(int i=0;i<(int)position.size();i++){position[i]=k[indice+i];}
}

void Corps::SubPosition(vector<double> k,double h,double indice){for(int i=0;i<(int)position.size();i++){position[i]=position[i]-k[indice+i]*h;}
}

void Corps::SetVitesse(vector<double> k,double indice){for(int i=0;i<(int)vitesse.size();i++){vitesse[i]=k[indice+i];}
}

