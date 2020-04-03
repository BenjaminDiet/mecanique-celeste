#include <iostream>
#include <cmath> 
#include <fstream> 
#include <vector> 
#include <string> 
#include <chrono>
#include <iomanip>

using namespace std;


#include "Utiles.h"
#include "Corps.h"
#include "Systeme.h"
#include "Constantes.h"
#include "Resolution.h"




int main(int argc, char *argv[]){

	double multi(stod(argv[3]));	// Nombre d'années
	n = atoi(argv[2]);		// Nombre de points de mesure
	string unit  = argv[4];
	string dataLink = argv[1];		// Fichier de vector initiaux
	int mS(3), mE(3);

	vector <double> sortiesActivees(4, true); // Affichage position, affichage aires, affichage energie, excentricité/périodes
	vector <double> sortiesRefusees(4, false); // Affichage position, affichage aires, affichage energie, excentricité/périodes
	vector <double> sortiesExcentricite={false, false, false, true}; 
	
	
	h = annee*multi/(n);

	if(string(argv[6]) == "relativiste") relativiste = true;

	if(unit == "UA") coeffPos = 1/(1.496e+11);

	if(string(argv[5]) == "Euler"){mS = 0 ; mE = 0;}
	if(string(argv[5]) == "Cromer"){mS = 1 ; mE = 1;}
	if(string(argv[5]) == "Verlet"){mS = 2 ; mE = 2;}
	if(string(argv[5]) == "RK4"){mS = 3 ; mE = 3;}
	if(string(argv[5]) == "all"){mS = 0 ; mE = 3;}



	
	Systeme systeme(dataLink,h);
	Systeme sys(dataLink,h);
	
	cout << endl <<"---- PARAMETRES ----"<< endl;
	cout << "\t" << "Programme " << dataLink << endl;
	cout << "\t" << systeme.getSize() << " objets." << endl;
	cout << "\t" << "Unite [" << unit << "]" << endl;
	cout << "\t" << multi << " ans." << endl;
	cout << "\t" << n << " points." << endl;
	cout << "\t" << h << " sec/pt ou " << (h/3600) << " h/pt" << endl;
	cout << "-- FIN PARAMETRES --"<< endl;
	cout << endl;


	int calcul = 0;

	
		
	
						
			if(calcul==0){

				for(idMethode = mS ; idMethode <= mE ; idMethode++){ // Choix méthodes
					nomMethode = methodes[idMethode];
					cout << "\tMethode " << nomMethode << endl;
					chrono::steady_clock::time_point s = chrono::steady_clock::now(); // Début chrono


					sys = resoudreSysteme(systeme, sortiesActivees); // Résoudre système dans le bon sens
				/* 
					// ALLER RETOUR
					h *= -1;
					sys = resoudreSysteme(sys, sortiesRefusees); // Résoudre dans l'autre sens	
					
					cout << "\t\tDistance apres aller-retour \t" << comparaisonAllerRetour(systeme.getPositions(), sys.getPositions())*coeffPos << " " << unit << endl;
				*/

		cout << "\t\tTemps méthode " << " \t\t\t" << chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - s).count()/1000 << "[ms]" << endl<<endl; // Fin chrono
	
				}
			}
			if(calcul==1){
				do{
					sys = resoudreSysteme(systeme, sortiesExcentricite); // Résoudre système dans le bon sens
					cout << norme(systeme[3].getVitesse()) << endl;
					systeme[3].multiplierVitesse(1.001);
				}while(sys[3].getExcentricite() > 0);
			}	

	




	return 0;
}
