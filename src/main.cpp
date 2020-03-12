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

	double annee(31536000); // 1 année
	double multi(stod(argv[3]));	// Nombre d'années
	double n(atoi(argv[2]));		// Nombre de points de mesure
	string unit  = argv[4];
	string dataLink = argv[1];		// Fichier de vector initiaux
	double coeffPos(1);
	int mS(3), mE(3);
	string methode[] = {"Euler", "Euler-Cromer", "Verlet", "RK4"};
	bool relativiste(false);
	

	if(string(argv[6]) == "relativiste") relativiste = true;
	if(string(argv[6]) == "classique") relativiste = false;

	double h = annee*multi/(n);
	if(unit == "UA") coeffPos = 1/(1.496e+11);

	if(string(argv[5]) == "Euler"){mS = 0 ; mE = 0;}
	if(string(argv[5]) == "Cromer"){mS = 1 ; mE = 1;}
	if(string(argv[5]) == "Verlet"){mS = 2 ; mE = 2;}
	if(string(argv[5]) == "RK4"){mS = 3 ; mE = 3;}
	if(string(argv[5]) == "all"){mS = 0 ; mE = 3;}





	Systeme systeme(dataLink);
	vector<vector <double>> coordInitiales = systeme.getPositions();


	cout << endl <<"-- PARAMETRES --"<< endl;
	cout << "\t" << "Programme " << dataLink << endl;
	cout << "\t" << systeme.getSize() << " objets." << endl;
	cout << "\t" << "Unite [" << unit << "]" << endl;
	cout << "\t" << multi << " ans." << endl;
	cout << "\t" << n << " points." << endl;
	cout << "\t" << h << " sec/pt ou " << (h/3600) << " h/pt" << endl;
	cout << endl;






	for(int methodeID = mS ; methodeID <= mE ; methodeID++){
		
		chrono::steady_clock::time_point s = chrono::steady_clock::now();
		
		Systeme sys(dataLink);
		
		
		cout << "\tMethode " << methode[methodeID] << endl;

		// Résoudre système dans le bon sens
		sys = resoudreSysteme(systeme, methodeID, n, h, coeffPos, relativiste);	
		// Résoudre système dans le temps négatif
		sys = resoudreSysteme(sys, methodeID, n, -h, coeffPos, relativiste);	

		// Pour comparaison des distances
		vector<vector <double>> coordFinales = sys.getPositions();
		double erreur(0);
		for(int i = 0 ; i < (int) coordFinales.size() ; i++){
			erreur += norme(distance(coordInitiales[i], coordFinales[i]));
		}
		erreur /= coordFinales.size();
		cout << "\t\tDistance apres aller-retour : " << erreur*coeffPos << " " << unit << endl;
		// Fin comparaison


		chrono::steady_clock::time_point f = chrono::steady_clock::now();		
		cout << "\t\tTemps aller-retour " << " = \t" << chrono::duration_cast<chrono::microseconds>(f - s).count()/1000 << "[ms]" << endl<<endl;
	}




	return 0;
}
