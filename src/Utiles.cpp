#include "Utiles.h"
#include "Corps.h"
#include "Systeme.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>

using namespace std;



void afficher(const vector<double>& v ){
	for(int i =0 ; i < (int) v.size() ; i++){cout << i << ") " << v[i] << endl;}
}

double norme(const vector<double>& v){
	double nor(0);
	for(int i = 0 ; i < (int) v.size() ; i++){
		nor += v[i]*v[i];
	}
	nor = sqrt(nor);
	return nor;
}


vector<double> normaliser(const vector<double>& v){
	vector<double> a;	
	a = multScalaire(1.0/norme(v), v);
	return a;
}

double produitScalaire(const vector<double>& a, const vector<double>& b){
	double s(0);	
	for(int i = 0 ; i < (int) a.size() ; i++){
		s+= a[i]*b[i];
	}
	return s;
}

double moyenne(const vector<double>& v){
	double moy;
	for(int i = 0 ; i < (int) v.size() ; i++){
		moy += v[i];
	}
	moy /= v.size();
	return moy;
}

vector<double> multScalaire(double m, const vector<double>& a){
	vector <double> v;
	for(int i = 0 ; i < (int) a.size() ; i++){
		v.push_back(a[i]*m);
	}
	return v;
}

vector<double> distance(const vector<double>& a, const vector<double>& b){
	vector<double> distance;
	if(a.size() == b.size()){
		for(int i = 0 ; i < (int) a.size() ; i++){
				distance.push_back(b[i]-a[i]);
		}
	}
	return distance;
}


vector<double> add(const vector<double>& a, const vector<double>& b){
	vector<double> v;
	if(a.size() == b.size()){
		for(int i = 0 ; i < (int) a.size() ; i++){
				v.push_back(b[i]+a[i]);
		}		
	}
	return v;
}
vector<double> oppose(const vector<double>& a){

	vector<double> oppose;
	oppose = multScalaire(-1, a);
	return oppose;
	
}

vector<double> ProdVec(const vector<double>& U,const vector<double>& V){
	
	vector<double> V3(U.size());
	
	V3[0]=U[1]*V[2]-U[2]*V[1];
	V3[1]=U[2]*V[0]-U[0]*V[2];
	V3[2]=U[0]*V[1]-U[1]*V[0];
	return V3;
}





Systeme resoudreSysteme(Systeme systeme, int id, int n, double h, double coeffPos, bool relativiste, vector <double> sorties){


	string methode;
	if(id == 0) methode = "Euler";
	if(id == 1) methode = "EulerCromer";
	if(id == 2) methode = "Verlet";
	if(id == 3) methode = "RK4";


	// STOCKAGE POUR ECRITURE
	vector<vector <vector <double>>> coordinates((int) systeme.getSize()); // Stockage des coordonnees pour ecriture
	vector<vector <double>> aires((int) systeme.getSize()); // Stockage des aires
	vector<double> airesInitiales; // Stockage des aires au temps zero
	vector<double> energieMeca; // Stockage des energies cinétiques


	systeme.calculerBarycentre();
	systeme.calculerAires(h);
	airesInitiales = multScalaire(pow(coeffPos,2),systeme.getAires());

	// CALCUL	
	double energie = 0;
	double energie0 = 0;
	double ecartAire = 0;



// === DEBUT RESOLUTION SYSTEME === //
	for(int k = 0 ; k < n ; k++){

		
		vector <double> aire = systeme.getAires();
		vector <vector <double>> pos = systeme.getPositions();
		
		// RESOLUTION
		if(id == 0) systeme.resoudreEuler(h, relativiste); // Résoud Euler à l'instant k*h pour tous les corps
		if(id == 1) systeme.resoudreEulerCromer(h, relativiste);
		if(id == 2) systeme.resoudreVerlet(h, relativiste);
		if(id == 3) systeme.resoudreRK4(h, relativiste);
		// FIN RESOLUTION
		// systeme.calculerBarycentre();
		// systeme.centrerBarycentre();
		systeme.calculerAires(h);


		energie = 0;
		ecartAire = 0;

		for(int i = 0 ; i < (int) pos.size(); i++){ // Tous les corps


			// CHANGEMENT UNITES
			vector <double> posT;
			for(int m = 0 ; m < (int) pos[i].size() ; m++){posT.push_back(pos[i][m]*coeffPos);}
			coordinates[i].push_back(posT); // ajoute les coordonnees a sortir
			// POSITIONS CONVERTIES

			
			// ECART AIRE PLANETE
			ecartAire = (aire[i]*coeffPos*coeffPos-airesInitiales[i])/airesInitiales[i]*100.0; // Energie devient erreur relative
			aires[i].push_back(ecartAire);
			// FIN ECART AIRE PLANETE


			// CALCUL ENERGIES
				// ENERGIE CIN
					energie += 0.5*systeme[i].getMasse()*pow(norme(systeme[i].getVitesse()),2);
				// ENERGIE MECA
			for(int j = i+1 ; j < (int) pos.size() ; j++) {
				vector <double> distanceVector = distance(systeme[i].getPosition(), systeme[j].getPosition());
				energie-=6.67e-11*systeme[i].getMasse()*systeme[j].getMasse()/norme(distanceVector);
			}
		}
		
		// ECART ENERGIE
		if(k==0)energie0 = energie; // Enregistre première énergie comme référence
		else{
			energie = (energie-energie0)/energie0*100.0; // Energie devient erreur relative
			energieMeca.push_back(energie);
		}
		// FIN ECART ENERGIE
		
	}

// === FIN RESOLUTION SYSTEME === //


// === TRAITEMENT FICHIER === //



	ofstream f;
	
	if(sorties[3]){
		for(int i = 1 ; i < (int)  coordinates.size() ; i++){ // Pour toutes les planètes
		
			if(systeme[i].getNature() == 1){ // Juste les planètes
				double diff = norme(coordinates[i][1]) - norme(coordinates[i][0]);
				double mini(0), maxi(0);
				int j = 1;
				if(diff >= 0){
					while((norme(coordinates[i][j+1]) > norme(coordinates[i][j])) and (j < ((int) coordinates[i].size()-2))) j++;
					maxi = j;
					while((norme(coordinates[i][j+1]) < norme(coordinates[i][j]))  and (j < ((int) coordinates[i].size()-2)))j++;
					mini = j;
				}
				else if(diff < 0){
					while((norme(coordinates[i][j+1]) < norme(coordinates[i][j])) and (j < ((int) coordinates[i].size()-2))) j++;
					mini = j;
					while((norme(coordinates[i][j+1]) > norme(coordinates[i][j])) and (j < ((int) coordinates[i].size()-2))) j++;
					maxi = j;		
				}
				double a = norme(coordinates[i][maxi]);
				double p = norme(coordinates[i][mini]);	
				double e = 1.0 - 2.0/(a/p + 1);
				double periode = abs(2.0*(maxi-mini)*h/(86400));				
				systeme[i].setPeriode(p);
				systeme[i].setExcentricite(e);
				
				f.open("../periodes/Periode"+methode+"_"+systeme[i].getLien(),fstream::app); // Récupère le lien
				f << periode<< endl;	
				f.close();	
				f.open("../eccentricites/Ecc"+methode+"_"+systeme[i].getLien(),fstream::app); // Récupère le lien
				f << e << endl;	
				f.close();	
			}
		}
	}	


	if(sorties[0]){

		
	
		// Ecriture des coordonnees dans des fichiers separes
			for(int i = 0 ; i < (int) coordinates.size() ; i++){ // Pour toutes les planètes
			f.open("../positions/"+methode+"_"+systeme[i].getLien(),fstream::app); // Récupère le lien
			for(int j = 0 ; j < (int) coordinates[i].size(); j++){ // tous les points
				f << setprecision(20) << (j*h)/(31557600) << "\t";
				for(int k = 0 ; k < (int) coordinates[i][j].size() ; k++){
					f << coordinates[i][j][k] << "\t";
				}
				f << endl;
			}

			f.close();
		}
		
		
		
		
	}
	
	
	if(sorties[1]){

		// Ecriture des aires dans un fichier global
		for(int i = 0 ; i < (int) aires.size() ; i++){ // Pour toutes les planètes 
			if(systeme[i].getNature() == 1){ // Juste les planètes
				f.open("../aires/"+methode+"_"+systeme[i].getLien(),fstream::app); // Récupère le lien
				for(int j = 0 ; j < (int) aires[i].size(); j++){
					f << setprecision(20) << (j*h)/(31557600) << "\t" << aires[i][j] << endl;
				}
				f.close();
			}
		}
	}


	if(sorties[2]){
		// Ecriture des énergies
		f.open("../energies/energieMecanique"+methode +".txt",fstream::app);
		for(int i = 0 ; i < (int) energieMeca.size() ; i++){
			f << setprecision(20) << (i*h)/(31557600) << "\t" << energieMeca[i] << endl;
		}
		f.close();
	}




	return systeme;
}





double comparaisonAllerRetour(vector<vector <double>> coordInitiales, vector<vector <double>> coordFinales){

	double erreur(0);
	for(int i = 0 ; i < (int) coordFinales.size() ; i++){
		erreur += norme(distance(coordInitiales[i], coordFinales[i]));
	}
	erreur /= coordFinales.size();
	return erreur;
}
