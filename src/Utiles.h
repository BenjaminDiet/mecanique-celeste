#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>

using namespace std;

class Systeme;

void afficher(const vector<double>& v ); // affiche un vector

double norme(const vector<double>& v); // norme d'un vector
double moyenne(const vector<double>& v); // moyenne d'un vector
vector<double> normaliser(const vector<double>& v);
vector<double> ProdVec(const vector<double>& U,const vector<double>& V);
double produitScalaire(const vector<double>& a, const vector<double>& b);
vector<double> distance(const vector<double>& a, const vector<double>& b);
vector<double> add(const vector<double>& a, const vector<double>& b);
vector<double> oppose(const vector<double>& a);
vector<double> multScalaire(double m, const vector<double>& a);


Systeme resoudreSysteme(Systeme systeme, int id, int n, double h, double coeffPos, bool relativiste);
