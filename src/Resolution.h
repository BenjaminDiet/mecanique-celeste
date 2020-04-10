#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include "parametre.h"
using namespace std;

class Systeme;



Systeme resoudreSysteme(Systeme systeme, parametre para);


void calculerEnergiesMecaniques(Systeme & sys, vector<vector <vector <double>>> & coord,vector<vector <vector <double>>> & vitesses,parametre para);
void calculerExcentricitesPeriodes(Systeme & sys, vector<vector <vector <double>>> & coord,parametre para);
void ecrirePositions(Systeme & sys, vector<vector <vector <double>>> & coord,parametre para);
void calculerAires(Systeme & sys, vector<vector <vector <double>>> & coord, vector<vector <vector <double>>> & vitesses, parametre para);

double comparaisonAllerRetour(vector<vector <double>> coordInitiales, vector<vector <double>> coordFinales, parametre para);




double calculerEnergiesMecaniques(Systeme & sys, int id, vector <vector <double>> & coord,vector <vector <double>> & vitesses);