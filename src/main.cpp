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



	string dataLink = string(argv[2]);
	int mS(3), mE(3);
	int calcul = 0;
	string unit = "";

	if(string(argv[1]) == "simuler"){ // Simulation du système
		calcul = 0;
		para.n = atoi(argv[3]);
		para.h = annee*stod(argv[4])/para.n;
		unit = string(argv[5]);
		if(unit == "UA") para.coeffPos = 1/(1.496e+11);
		
		if(string(argv[6]) == "Euler"){mS = 0 ; mE = 0;para.nomMethode="Euler";para.idMethode=0;}
		if(string(argv[6]) == "Cromer"){mS = 1 ; mE = 1;para.nomMethode="Cromer";para.idMethode=1;}
		if(string(argv[6]) == "Verlet"){mS = 2 ; mE = 2;para.nomMethode="Verlet";para.idMethode=2;}
		if(string(argv[6]) == "RK4"){mS = 3 ; mE = 3;para.nomMethode="RK4";para.idMethode=3;}
		if(string(argv[6]) == "all"){mS = 0 ; mE = 3;para.nomMethode="all";para.idMethode=0;}
	

		if(string(argv[7]) == "relativiste") para.relativiste = true;

		if(string(argv[8]) == "reversible") para.reversible = true;

		if(atoi(argv[9])==1)  para.sorties[0] = true;
		if(atoi(argv[10])==1) para.sorties[1] = true;
		if(atoi(argv[11])==1) para.sorties[2] = true;
		if(atoi(argv[12])==1) para.sorties[3] = true;
	
	}else if(string(argv[1]) == "liberation"){ // Vitesse de libération
		calcul = 1;
		para.idPlanete = atoi(argv[3]);
		para.precision = pow(10, -atoi(argv[4]));
		cout << para.idPlanete << endl;
		cout << para.precision << endl;
	}




	
	Systeme systeme(dataLink, para.h);
	Systeme sys(dataLink, para.h);

	if(calcul == 0){
		cout << endl << "|----- PARAMETRES CALCUL -----|"<< endl;
		cout << "|" << "Programme \"" << dataLink << "\"" << endl;
		cout << "|" << "Méthode : " << para.nomMethode  << endl;
		if(para.relativiste) cout << "|" << "Approche classique." << endl;
		else cout << "|" << "Approche relativiste." << endl;
		cout << "|" << "Réversibilité : " << para.reversible  << endl;
		cout << "|" << systeme.getSize() << " objets." << endl;
		cout << "|" << "Unite [" << string(argv[5]) << "]" << endl;
		cout << "|" << stod(argv[4]) << " ans." << endl;
		cout << "|" << para.n << " points." << endl;
		cout << "|" << para.h << " sec/pt ou " << (para.h/3600) << " h/pt" << endl;
		cout << "|------ PARAMETRES AFFICHAGE ------|"<< endl;
		cout << "|Positions : " << para.sorties[0] << endl;
		cout << "|Aires : " << para.sorties[1] << endl;
		cout << "|Energie Meca : " << para.sorties[2] << endl;
		cout << "|Excentricités/périodes : " << para.sorties[3] << endl;
		cout << "|--------------------------------|"<< endl;
		cout << endl;
	}
	if(calcul == 1){
		cout << endl << "|----- PARAMETRES CALCUL -----|"<< endl;
		cout << "|" << "Programme \"" << dataLink << "\"" << endl;
		cout << "|" << systeme.getSize() << " objets." << endl;
		cout << "|" << "Planète : " << systeme[para.idPlanete].getNom() << "." << endl;
		cout << "|" << "Précision : " << para.precision  << endl;
		cout << endl;
	}


			if(calcul==1){

				vector <double> a, b, m;
				double fA(-1), fB(1), fM(0); // energies mécanique vitesses a, b, et m
				vector <vector <double>> pos, vit;

				a = sys[(int) para.idPlanete].getVitesse();// CALCUL ENERGIE MECA en A
					pos = sys.getPositions();
					vit = sys.getVitesses();
					fA = calculerEnergiesMecaniques(sys, para.idPlanete, pos, vit);
				// FIN CALCUL				
				b = multScalaire(10,sys[para.idPlanete].getVitesse());
				sys[para.idPlanete].SetVitesse(b,0);
				// CALCUL ENERGIE MECA en B
					pos = sys.getPositions();
					vit = sys.getVitesses();
					fB = calculerEnergiesMecaniques(sys, para.idPlanete, pos, vit);
				// FIN CALCUL

			do{	

				m = multScalaire(0.5,(a+b));
				sys[para.idPlanete].SetVitesse(m,0);
				// CALCUL ENERGIE MECA en B
					pos = sys.getPositions();
					vit = sys.getVitesses();
					fM = calculerEnergiesMecaniques(sys, para.idPlanete, pos, vit);
				// FIN CALCUL



				if(fM < 0){a = m; fA = fM;}
				if(fM > 0){b = m; fB = fM;}
			}while(norme(a-b) > para.precision);



			cout << "Vitesse de libération : " << setprecision(20) << norme(sys[para.idPlanete].getVitesse()) << endl;


			}



			if(calcul==0){

				for(para.idMethode = mS ; para.idMethode <= mE ; para.idMethode++){ // Choix méthodes
					para.nomMethode = para.methodes[para.idMethode];
					cout << "\tMethode " << para.nomMethode << endl;
					chrono::steady_clock::time_point s = chrono::steady_clock::now(); // Début chrono


					sys = resoudreSysteme(systeme, para); // Résoudre système dans le bon sens
				 if(para.reversible){
					para.h *= -1;
					para.sorties[0] = false;
					para.sorties[1] = false;
					para.sorties[2] = false;
					para.sorties[3] = false;
					sys = resoudreSysteme(sys, para); // Résoudre dans l'autre sens	
					
					cout << "\t\tDistance apres aller-retour \t" << comparaisonAllerRetour(systeme.getPositions(), sys.getPositions(),para)*para.coeffPos << " " << unit << endl;
				}

				cout << "\t\tTemps méthode " << " \t\t\t" << chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - s).count()/1000 << "[ms]" << endl<<endl; // Fin chrono
	
				}
			}

				

	




	return 0;
}
