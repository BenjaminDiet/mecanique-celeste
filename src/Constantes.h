#include <iostream>
#include <string>


// CONSTANTES PHYSIQUES
extern double G;
extern double c;
extern double annee;

// PARAMETRES DE SIMULATION
extern double coeffPos;	// Changement d'unité en sortie
extern double h;	// dt
extern bool relativiste;
extern int n; // Nombre de points


// METHODES DE RESOLUTION
extern std::string nomMethode;
extern int idMethode;
extern std::string methodes[];