#ifndef DEF_PARAMETRE
#define DEF_PARAMETRE
 
#include <string>
 
struct parametre{
double coeffPos = 1;
double h = 0;
bool relativiste = false;
int n = 1;
std::string nomMethode = "";
int idMethode = 0;
std::string methodes[4] =  {"Euler", "EulerCromer", "Verlet", "RK4"};	
};
 
#endif
