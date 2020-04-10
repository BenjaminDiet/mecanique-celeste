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
#include "parametre.h"




int main(int argc, char *argv[]){
	parametre para;
	double multi(stod(argv[3]));	// Nombre d'années
	para.n = atoi(argv[2]);		// Nombre de points de mesure
	string unit  = argv[4];
	string dataLink = argv[1];		// Fichier de vector initiaux
	int mS(3), mE(3);

	vector <double> sortiesActivees(4, true); // Affichage position, affichage aires, affichage energie, excentricité/périodes
	vector <double> sortiesRefusees(4, false); // Affichage position, affichage aires, affichage energie, excentricité/périodes
	vector <double> sortiesExcentricite={false, false, false, true}; 
	
	
	para.h = annee*multi/(para.n);

	if(string(argv[6]) == "relativiste") para.relativiste = true;

	if(unit == "UA") para.coeffPos = 1/(1.496e+11);

	if(string(argv[5]) == "Euler"){mS = 0 ; mE = 0;para.nomMethode="Euler";para.idMethode=0;}
	if(string(argv[5]) == "Cromer"){mS = 1 ; mE = 1;para.nomMethode="Cromer";para.idMethode=1;}
	if(string(argv[5]) == "Verlet"){mS = 2 ; mE = 2;para.nomMethode="Verlet";para.idMethode=2;}
	if(string(argv[5]) == "RK4"){mS = 3 ; mE = 3;para.nomMethode="RK4";para.idMethode=3;}
	if(string(argv[5]) == "all"){mS = 0 ; mE = 3;para.nomMethode="all";para.idMethode=0;}
	


	
	Systeme systeme(dataLink,para.h);
	Systeme sys(dataLink,para.h);
	
	cout << endl <<"---- PARAMETRES ----"<< endl;
	cout << "\t" << "Programme " << dataLink << endl;
	cout << "\t" << systeme.getSize() << " objets." << endl;
	cout << "\t" << "Unite [" << unit << "]" << endl;
	cout << "\t" << multi << " ans." << endl;
	cout << "\t" << para.n << " points." << endl;
	cout << "\t" << para.h << " sec/pt ou " << (para.h/3600) << " h/pt" << endl;
	cout << "-- FIN PARAMETRES --"<< endl;
	cout << endl;


	int calcul = 0;

	double precision = 1e-6;

	int id = 3;

			if(calcul==1){

				vector <double> a, b, m;
				double fA(-1), fB(1), fM(0); // energies mécanique vitesses a, b, et m
				vector <vector <double>> pos, vit;

				a = sys[id].getVitesse();// CALCUL ENERGIE MECA en A
					pos = sys.getPositions();
					vit = sys.getVitesses();
					fA = calculerEnergiesMecaniques(sys, id, pos, vit);
				// FIN CALCUL				
				b = multScalaire(10,sys[id].getVitesse());
				sys[id].SetVitesse(b,0);
				// CALCUL ENERGIE MECA en B
					pos = sys.getPositions();
					vit = sys.getVitesses();
					fB = calculerEnergiesMecaniques(sys, id, pos, vit);
				// FIN CALCUL

			do{	

				m = multScalaire(0.5,(a+b));
				cout << setprecision(10) << norme(m) << endl;
				sys[id].SetVitesse(m,0);
				// CALCUL ENERGIE MECA en B
					pos = sys.getPositions();
					vit = sys.getVitesses();
					fM = calculerEnergiesMecaniques(sys, id, pos, vit);
				// FIN CALCUL



				if(fM < 0){a = m; fA = fM;}
				if(fM > 0){b = m; fB = fM;}
			}while(norme(a-b) > precision);



			cout << "Vitesse : " << norme(sys[id].getVitesse()) << "\t";


			}



			if(calcul==0){

				for(para.idMethode = mS ; para.idMethode <= mE ; para.idMethode++){ // Choix méthodes
					para.nomMethode = para.methodes[para.idMethode];
					cout << "\tMethode " << para.nomMethode << endl;
					chrono::steady_clock::time_point s = chrono::steady_clock::now(); // Début chrono


					sys = resoudreSysteme(systeme, sortiesActivees,para); // Résoudre système dans le bon sens
				 /*
					// ALLER RETOUR
					para.h *= -1;
					sys = resoudreSysteme(sys, sortiesRefusees,para); // Résoudre dans l'autre sens	
					
					cout << "\t\tDistance apres aller-retour \t" << comparaisonAllerRetour(systeme.getPositions(), sys.getPositions(),para)*para.coeffPos << " " << unit << endl;
				*/

				cout << "\t\tTemps méthode " << " \t\t\t" << chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - s).count()/1000 << "[ms]" << endl<<endl; // Fin chrono
	
				}
			}

				

	




	return 0;
}
