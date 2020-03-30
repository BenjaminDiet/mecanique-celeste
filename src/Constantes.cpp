#include "Constantes.h"
#include <string>

using namespace std;


double G = 6.6743015e-11;
double c = 299792458;
double annee = (31557600); // 1 année en secondes


double coeffPos = 1;
double h = 0;
bool relativiste = false;
int n = 1;


string nomMethode = "";
int idMethode = 0;
string methodes[] =  {"Euler", "EulerCromer", "Verlet", "RK4"};