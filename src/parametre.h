#ifndef DEF_PARAMETRE
#define DEF_PARAMETRE
 
#include <string>
#include <vector>
 
struct parametre{
	// Simulation
	bool reversible = false;
	double coeffPos = 1;
	double h = 1;
	bool relativiste = false;
	int n = 1;
	std::string nomMethode = "";
	int idMethode = 0;
	std::vector <bool> sorties = {false, false, false, false};
	std::string methodes[4] =  {"Euler", "EulerCromer", "Verlet", "RK4"};	

	// Vitesse de libération
	int idPlanete = 1;
	double precision = 1;
};
 
#endif
