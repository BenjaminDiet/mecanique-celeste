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



Systeme resoudreSysteme(Systeme systeme, vector <double> sorties,parametre para){


	// STOCKAGE POUR ECRITURE
	vector<vector <vector <double>>> coordinates((int) systeme.getSize()); // Stockage des coordonnees pour ecriture
	vector<vector <vector <double>>> velocities((int) systeme.getSize()); // Stockage des vitesses pour ecriture


	

// === DEBUT RESOLUTION SYSTEME === //
	for(int k = 0 ; k <= para.n ; k++){		
		// RESOLUTION
		if(para.idMethode == 0) systeme.resoudreEuler(para.h, para.relativiste); // Résoud Euler à l'instant k*h pour tous les corps
		if(para.idMethode == 1) systeme.resoudreEulerCromer(para.h, para.relativiste);
		if(para.idMethode == 2) systeme.resoudreVerlet(para.h, para.relativiste);
		if(para.idMethode == 3) systeme.resoudreRK4(para.h, para.relativiste);
		// FIN RESOLUTION



		

		for(int i = 0 ; i < systeme.getSize(); i++){ // Tous les corps

			// ENREGISTREMENT POSITIONS ET VITESSES ET AIRES
			coordinates[i].push_back(systeme[i].getPosition()); // enregistre les coordonnees a sortir
			velocities[i].push_back(systeme[i].getVitesse()); // enregistre les vitesses a sortir
			// POSITIONS CONVERTIES ET ENREGISTREES	
		}


	}

// === FIN RESOLUTION SYSTEME === //

// === TRAITEMENT FICHIER === //



	
	if(sorties[3]){	calculerExcentricitesPeriodes(systeme, coordinates,para);}	


	// ECRITURE POSITIONS
	if(sorties[0]) ecrirePositions(systeme, coordinates, para);
	

	// ECRITURE AIRES
	if(sorties[1]) ecrireAires(systeme, coordinates, velocities,para);

	// ECRITURE ENERGIE MECA
	if(sorties[2]) calculerEnergiesMecaniques(systeme, coordinates, velocities,para);


	return systeme;
}







void ecrireAires(Systeme & sys, vector<vector <vector <double>>> & coord, vector<vector <vector <double>>> & vitesses,parametre para){
	vector<vector <double>> airesTotales;


	// CALCUL AIRES
	for(int t = 0 ; t < (int) coord[0].size() ; t++){	// Chaque instant t 
		vector <double> airesT;
		for(size_t i = 1 ; i < coord.size() ; i++){	// Chaque planète i sauf Soleil
			double aire = 0;
			aire = norme(ProdVec(distance(coord[i][t], coord[0][t]), vitesses[i][t]))*para.h/2.0;
			airesT.push_back(aire);
		}
		airesTotales.push_back(airesT);
	// FIN ENERGIES MECANIQUES
	}



	ofstream f;
	for(int i = 1 ; i < sys.getSize() ; i++){ // Pour toutes les planètes
			if(sys[i].getNature()==1){
			f.open("../aires/"+para.nomMethode+"_"+sys[i].getLien(),fstream::app); // Récupère le lien
		for(size_t t = 1 ; t < airesTotales.size(); t++){ // A chaque instant
			f << setprecision(20) << (t*para.h)/(annee) << "\t" <<(airesTotales[t][i]-airesTotales[0][i])/airesTotales[0][i]*100.0 << endl;
		}
		f.close();
	}
	}

}





void ecrirePositions(Systeme & sys, vector<vector <vector <double>>> & coord,parametre para){

	ofstream f;
 
	for(size_t i = 0 ; i < coord.size() ; i++){ // Pour toutes les planètes
			// OUVRE LE FICHIER
			f.open("../positions/"+para.nomMethode+"_"+sys[i].getLien(),fstream::app);

			for(size_t t = 0 ; t < coord[i].size(); t++){ // CHAQUE TEMPS T
				f << setprecision(20) << (t*para.h)/(annee) << "\t";
				for(int k = 0 ; k < (int) coord[i][t].size() ; k++){ // CHAQUE COORDONNEE
					f << coord[i][t][k]*para.coeffPos << "\t";
				}
				f << endl;
			}

			f.close();
		}
		
		
}







void calculerExcentricitesPeriodes(Systeme & sys, vector<vector <vector <double>>> & coord,parametre para){

	ofstream f;
	cout << "\t\tPlanète\tExcentricités \tPériodes" << endl;

	// CHAQUE PLANETE
	for(size_t i = 1 ; i < coord.size() ; i++){
			
		if(sys[i].getNature() == 1){ // JUSTE PLANETES
			double diff = norme(coord[i][1]) - norme(coord[i][0]);
			double mini(0), maxi(0);
			int j = 1;

			double n1 = norme(distance(coord[0][j+1],coord[i][j+1]));
			double n2 = norme(distance(coord[0][j+1],coord[i][j]));

			if(diff >= 0){
				while((n1 > n2) and (j < ((int) coord[i].size()-2))){
					j++;
					n1 = norme(distance(coord[0][j+1],coord[i][j+1]));
					n2 = norme(distance(coord[0][j+1],coord[i][j]));
				}
				maxi = j;
				while((n1 < n2)  and (j < ((int) coord[i].size()-2))){
					j++;
					n1 = norme(distance(coord[0][j+1],coord[i][j+1]));
					n2 = norme(distance(coord[0][j+1],coord[i][j]));
				}
				mini = j;
			}
			else if(diff < 0){
				while((n1 < n2) and (j < ((int) coord[i].size()-2))) {
					j++;
					n1 = norme(distance(coord[0][j+1],coord[i][j+1]));
					n2 = norme(distance(coord[0][j+1],coord[i][j]));
				}
				mini = j;
				while((n1 > n2) and (j < ((int) coord[i].size()-2))) {
					j++;
					n1 = norme(distance(coord[0][j+1],coord[i][j+1]));
					n2 = norme(distance(coord[0][j+1],coord[i][j]));
				}
				maxi = j;		
			}


			double a = norme(distance(coord[0][maxi],coord[i][maxi]));
			double p = norme(distance(coord[0][mini],coord[i][mini]));	
			double e = 1.0 - 2.0/(a/p + 1);
			double periode = abs(2.0*(maxi-mini)*para.h/(86400));	

			sys[i].setPeriode(periode);
			sys[i].setExcentricite(e);

			// Ecriture
			cout << "\t\t" << sys[i].getNom() << "\t" << e << "\t" << periode << endl;
			f.open("../periodes/Periode"+para.nomMethode+"_"+sys[i].getLien(),fstream::app); // Récupère le lien
			f << periode<< endl;	
			f.close();	
			f.open("../excentricites/Ecc"+para.nomMethode+"_"+sys[i].getLien(),fstream::app); // Récupère le lien
			f << e << endl;	
			f.close();	
		}
	}


}




void calculerEnergiesMecaniques(Systeme & sys, vector<vector <vector <double>>> & coord, vector<vector <vector <double>>> & vitesses,parametre para){
	
	vector<double> energiesMecaniques;
	double energieMeca(0);


	// CALCUL ENERGIES MECANIQUES
	for(int t = 0 ; t < (int) coord[0].size() ; t++){	// Chaque instant t 
		for(size_t i = 0 ; i < coord.size() ; i++){	// Chaque planète i

			// ENERGIE CINETIQUE
			energieMeca += 0.5*sys[i].getMasse() * norme(vitesses[i][t]) * norme(vitesses[i][t]);
			// FIN ENERGIE CINETIQUE

			// ENERGIE MECA
			for(size_t j=i+1 ; j < coord.size() ; j++){	// Chaque planète j
				double distanceV = norme(distance(coord[i][t], coord[j][t]));
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
	f.open("../energies/energieMecanique"+para.nomMethode+".txt",fstream::app);

	// TRANSFORMATION EN ERREUR RELATIVE
	for(size_t i=1 ; i < energiesMecaniques.size() ; i++){
		f << setprecision(20) << (i*para.h)/(31557600) << "\t" << (energiesMecaniques[i] - energiesMecaniques[0])/energiesMecaniques[0]*100.0 << endl;
	}
	f.close();

	// Pour affichage
	cout << "\t\t" << "Energie Méca finale (%) \t" << abs((energiesMecaniques[energiesMecaniques.size()-1] - energiesMecaniques[0])/energiesMecaniques[0]*100.0) << endl;

}









double comparaisonAllerRetour(vector<vector <double>> coordInitiales, vector<vector <double>> coordFinales,parametre para){

	double erreur(0);
	for(size_t i = 0 ; i < coordFinales.size() ; i++){
		erreur += norme(distance(coordInitiales[i], coordFinales[i]));
	}
	erreur /= coordFinales.size();
	erreur *= para.coeffPos;
	return erreur;
}



// Renvoie l'énergie mécanique au dernier point POUR VITESSE LIBERATION
double calculerEnergiesMecaniques(Systeme & sys, int id, vector <vector <double>> & coord,vector <vector <double>> & vitesses){
	
	double energieMeca(0);
 

			// ENERGIE CINETIQUE
			energieMeca += 0.5*sys[id].getMasse() * norme(vitesses[id]) * norme(vitesses[id]);
			// FIN ENERGIE CINETIQUE

			// ENERGIE MECA
			for(size_t j=0 ; j < coord.size() ; j++){	// Chaque planète j
				if((int) j!=id){
					double distanceV = norme(distance(coord[id], coord[j]));
					energieMeca -= G * sys[id].getMasse() * sys[j].getMasse() / distanceV;
				}

			}
	


	return energieMeca;
}

