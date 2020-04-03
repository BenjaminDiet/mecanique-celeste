#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>

using namespace std;

class Systeme;



Systeme resoudreSysteme(Systeme systeme, vector <double> sorties);


void calculerEnergiesMecaniques(Systeme & sys, vector<vector <vector <double>>> & coord,vector<vector <vector <double>>> & vitesses);
void calculerExcentricitesPeriodes(Systeme & sys, vector<vector <vector <double>>> & coord);
void ecrirePositions(Systeme & sys, vector<vector <vector <double>>> & coord);
void ecrireAires(Systeme & sys, vector<vector <vector <double>>> & coord, vector<vector <vector <double>>> & vitesses);

double comparaisonAllerRetour(vector<vector <double>> coordInitiales, vector<vector <double>> coordFinales);
