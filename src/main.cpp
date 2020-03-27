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




int main(int argc, char *argv[]){

	double annee(31557600); // 1 année
	double multi(stod(argv[3]));	// Nombre d'années
	double n(atoi(argv[2]));		// Nombre de points de mesure
	string unit  = argv[4];
	string dataLink = argv[1];		// Fichier de vector initiaux
	double coeffPos(1);
	int mS(3), mE(3);
	string methode[] = {"Euler", "Euler-Cromer", "Verlet", "RK4"};
	bool relativiste(false);
	vector <double> sortiesActivees(4, true); // Affichage position, affichage aires, affichage energie, excentricité/périodes
	vector <double> sortiesRefusees(4, false); // Affichage position, affichage aires, affichage energie, excentricité/périodes
	vector <double> sortiesExcentricite={false, false, false, true}; 
	
	
	double h = annee*multi/(n);

	if(string(argv[6]) == "relativiste") relativiste = true;
	if(string(argv[6]) == "classique") relativiste = false;

	if(unit == "UA") coeffPos = 1/(1.496e+11);

	if(string(argv[5]) == "Euler"){mS = 0 ; mE = 0;}
	if(string(argv[5]) == "Cromer"){mS = 1 ; mE = 1;}
	if(string(argv[5]) == "Verlet"){mS = 2 ; mE = 2;}
	if(string(argv[5]) == "RK4"){mS = 3 ; mE = 3;}
	if(string(argv[5]) == "all"){mS = 0 ; mE = 3;}



	
	Systeme systeme(dataLink);
	Systeme sys(dataLink);
	
	cout << endl <<"-- PARAMETRES --"<< endl;
	cout << "\t" << "Programme " << dataLink << endl;
	cout << "\t" << systeme.getSize() << " objets." << endl;
	cout << "\t" << "Unite [" << unit << "]" << endl;
	cout << "\t" << multi << " ans." << endl;
	cout << "\t" << n << " points." << endl;
	cout << "\t" << h << " sec/pt ou " << (h/3600) << " h/pt" << endl;
	cout << endl;


	int calcul = 0;

	
		
		
	
						
			if(calcul==0){

				for(int methodeID = mS ; methodeID <= mE ; methodeID++){ // Choix méthodes
					cout << "\tMethode " << methode[methodeID] << endl;
					chrono::steady_clock::time_point s = chrono::steady_clock::now(); // Début chrono


					sys = resoudreSysteme(systeme, methodeID, n, h, coeffPos, relativiste, sortiesActivees); // Résoudre système dans le bon sens
					sys = resoudreSysteme(sys, methodeID, n, -h, coeffPos, relativiste, sortiesRefusees); // Résoudre dans l'autre sens	
					
					cout << "\t\tDistance apres aller-retour : " << comparaisonAllerRetour(systeme.getPositions(), sys.getPositions())*coeffPos << " " << unit << endl;
					

		cout << "\t\tTemps méthode " << " = \t" << chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - s).count()/1000 << "[ms]" << endl<<endl; // Fin chrono
	
				}
			}
			if(calcul==1){
				do{
					sys = resoudreSysteme(systeme, mS, n, h, coeffPos, relativiste, sortiesExcentricite); // Résoudre système dans le bon sens
					cout << norme(systeme[1].getVitesse()) << endl;
					systeme[1].multiplierVitesse(1.001);
				}while(sys[1].getExcentricite() > 0);
			}	

	




	return 0;
}
