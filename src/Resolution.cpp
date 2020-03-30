#include "Utiles.h"
#include "Corps.h"
#include "Systeme.h"
#include "Constantes.h"
#include "Resolution.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>

using namespace std;



Systeme resoudreSysteme(Systeme systeme, vector <double> sorties){



	// STOCKAGE POUR ECRITURE
	vector<vector <vector <double>>> coordinates((int) systeme.getSize()); // Stockage des coordonnees pour ecriture
	vector<vector <vector <double>>> velocities((int) systeme.getSize()); // Stockage des vitesses pour ecriture
	vector<vector <double>> aires((int) systeme.getSize()); // Stockage des aires
	vector<double> energieMeca; // Stockage des energies mécaniques


	systeme.calculerBarycentre();

	

// === DEBUT RESOLUTION SYSTEME === //
	for(int k = 0 ; k <= n ; k++){		
		// RESOLUTION
		if(idMethode == 0) systeme.resoudreEuler(h, relativiste); // Résoud Euler à l'instant k*h pour tous les corps
		if(idMethode == 1) systeme.resoudreEulerCromer(h, relativiste);
		if(idMethode == 2) systeme.resoudreVerlet(h, relativiste);
		if(idMethode == 3) systeme.resoudreRK4(h, relativiste);
		// FIN RESOLUTION


		systeme.calculerAires();	// Met à jour les aires du systèmes

		

		for(int i = 0 ; i < systeme.getSize(); i++){ // Tous les corps

			// ENREGISTREMENT POSITIONS ET VITESSES ET AIRES
			coordinates[i].push_back(multScalaire(coeffPos, systeme[i].getPosition())); // enregistre les coordonnees a sortir
			velocities[i].push_back(multScalaire(coeffPos, systeme[i].getVitesse())); // enregistre les vitesses a sortir
			if(sorties[1]) if(systeme[i].getNature()==1)aires[i].push_back(systeme.getAires()[i]);
			// POSITIONS CONVERTIES ET ENREGISTREES	
		}


	}

// === FIN RESOLUTION SYSTEME === //

// === TRAITEMENT FICHIER === //



	
	if(sorties[3]){	calculerExcentricitesPeriodes(systeme, coordinates);}	


	// ECRITURE POSITIONS
	if(sorties[0]) ecrirePositions(systeme, coordinates);
	

	// ECRITURE AIRES
	if(sorties[1]) ecrireAires(systeme, aires);

	// ECRITURE ENERGIE MECA
	if(sorties[2]) calculerEnergiesMecaniques(systeme, coordinates, velocities);


	return systeme;
}







void ecrireAires(Systeme & sys, vector<vector <double>> & aires){
	ofstream f;
		for(size_t i = 0 ; i < aires.size() ; i++){ // Pour toutes les planètes 
			if(sys[i].getNature() == 1){ // Juste les planètes
				f.open("../aires/"+nomMethode+"_"+sys[i].getLien(),fstream::app); // Récupère le lien
				for(size_t t = 1 ; t < aires[i].size(); t++){
					f << setprecision(20) << (t*h)/(annee) << "\t" << (aires[i][t]-aires[i][0])/aires[i][0] << endl;
				}
				f.close();
			}
		}

}





void ecrirePositions(Systeme & sys, vector<vector <vector <double>>> & coord){

	ofstream f;
 
	for(size_t i = 0 ; i < coord.size() ; i++){ // Pour toutes les planètes
			// OUVRE LE FICHIER
			f.open("../positions/"+nomMethode+"_"+sys[i].getLien(),fstream::app);

			for(size_t t = 0 ; t < coord[i].size(); t++){ // CHAQUE TEMPS T
				f << setprecision(20) << (t*h)/(annee) << "\t";
				for(int k = 0 ; k < (int) coord[i][t].size() ; k++){ // CHAQUE COORDONNEE
					f << coord[i][t][k] << "\t";
				}
				f << endl;
			}

			f.close();
		}
		
		
}







void calculerExcentricitesPeriodes(Systeme & sys, vector<vector <vector <double>>> & coord){

	ofstream f;
	cout << "\t\tPlanète\tExcentricités \tPériodes" << endl;

	// CHAQUE PLANETE
	for(size_t i = 1 ; i < coord.size() ; i++){
			
		if(sys[i].getNature() == 1){ // JUSTE PLANETES
			double diff = norme(coord[i][1]) - norme(coord[i][0]);
			double mini(0), maxi(0);
			int j = 1;
			if(diff >= 0){
				while((norme(coord[i][j+1]) > norme(coord[i][j])) and (j < ((int) coord[i].size()-2))) j++;
				maxi = j;
				while((norme(coord[i][j+1]) < norme(coord[i][j]))  and (j < ((int) coord[i].size()-2)))j++;
				mini = j;
			}
			else if(diff < 0){
				while((norme(coord[i][j+1]) < norme(coord[i][j])) and (j < ((int) coord[i].size()-2))) j++;
				mini = j;
				while((norme(coord[i][j+1]) > norme(coord[i][j])) and (j < ((int) coord[i].size()-2))) j++;
				maxi = j;		
			}
			double a = norme(coord[i][maxi]);
			double p = norme(coord[i][mini]);	
			double e = 1.0 - 2.0/(a/p + 1);
			double periode = abs(2.0*(maxi-mini)*h/(86400));				
			sys[i].setPeriode(periode);
			sys[i].setExcentricite(e);
			cout << "\t\t" << sys[i].getNom() << "\t" << e << "\t" << periode << endl;
			f.open("../periodes/Periode"+nomMethode+"_"+sys[i].getLien(),fstream::app); // Récupère le lien
			f << periode<< endl;	
			f.close();	
			f.open("../excentricites/Ecc"+nomMethode+"_"+sys[i].getLien(),fstream::app); // Récupère le lien
			f << e << endl;	
			f.close();	
		}
	}


}




void calculerEnergiesMecaniques(Systeme & sys, vector<vector <vector <double>>> & coord, vector<vector <vector <double>>> & vitesses){
	
	vector<double> energiesMecaniques;
	double energieMeca(0);


	// CALCUL ENERGIES MECANIQUES
	for(int t = 0 ; t < (int) coord[0].size() ; t++){	// Chaque instant t 
		for(size_t i = 0 ; i < coord.size() ; i++){	// Chaque planète i

			// ENERGIE CINETIQUE
			energieMeca += 0.5*sys[i].getMasse() * norme(vitesses[i][t]) * norme(vitesses[i][t]) / (coeffPos*coeffPos);
			// FIN ENERGIE CINETIQUE

			// ENERGIE MECA
			for(size_t j=i+1 ; j < coord.size() ; j++){	// Chaque planète j
				double distanceV = norme(distance(coord[i][t], coord[j][t]))/coeffPos;
				energieMeca -= G * sys[i].getMasse() * sys[j].getMasse() / distanceV;

			}
			// FIN ENERGIE MECA
		}
		energiesMecaniques.push_back(energieMeca);
		energieMeca = 0;
	}
	// FIN ENERGIES MECANIQUES


	// ECRITURE
	ofstream f;
	f.open("../energies/energieMecanique"+nomMethode+".txt",fstream::app);

	// TRANSFORMATION EN ERREUR RELATIVE
	for(size_t i=1 ; i < energiesMecaniques.size() ; i++){
		f << setprecision(20) << (i*h)/(31557600) << "\t" << (energiesMecaniques[i] - energiesMecaniques[0])/energiesMecaniques[0]*100.0 << endl;
	}
	f.close();

	// Pour affichage
	cout << "\t\t" << "Energie Méca finale (%) \t" << abs((energiesMecaniques[energiesMecaniques.size()-1] - energiesMecaniques[0])/energiesMecaniques[0]*100.0) << endl;

}









double comparaisonAllerRetour(vector<vector <double>> coordInitiales, vector<vector <double>> coordFinales){

	double erreur(0);
	for(size_t i = 0 ; i < coordFinales.size() ; i++){
		erreur += norme(distance(coordInitiales[i], coordFinales[i]));
	}
	erreur /= coordFinales.size();
	return erreur;
}
