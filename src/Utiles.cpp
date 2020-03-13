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
	double nor;
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

	for(int k = 0 ; k < n ; k++){

		
		vector <double> aire = systeme.getAires();

		if(id == 0) systeme.resoudreEuler(h, relativiste); // Résoud Euler à l'instant k*h pour tous les corps
		if(id == 1) systeme.resoudreEulerCromer(h, relativiste);
		if(id == 2) systeme.resoudreVerlet(h, relativiste);
		if(id == 3) systeme.resoudreRK4(h, relativiste);
		systeme.calculerBarycentre();
		systeme.calculerAires(h);

		vector <vector <double>> pos = systeme.getPositions();

		energie = 0;
		ecartAire = 0;

		for(int i = 0 ; i < (int) pos.size(); i++){ // Tous les corps


			//Changement unité
			vector <double> posT;

			for(int k = 0 ; k < (int) pos[i].size() ; k++){posT.push_back(pos[i][k]*coeffPos);}
			coordinates[i].push_back(posT); // ajoute les coordonnees a sortir


			
			ecartAire = (aire[i]*coeffPos*coeffPos-airesInitiales[i])/airesInitiales[i]*100.0; // Energie devient erreur relative
			aires[i].push_back(ecartAire);


			// Calcul énergies
			// Cinétique
			energie += 0.5*systeme.getObjet(i).getMasse()*pow(norme(systeme.getObjet(i).getVitesse()),2);
			// Mécanique
			for(int j = i+1 ; j < (int) pos.size() ; j++) {
					vector <double> distanceVector(3);
					distanceVector=distance(systeme.getObjet(i).getPosition(), systeme.getObjet(j).getPosition());
					// Energie mécanique
					energie-=6.67e-11*systeme.getObjet(i).getMasse()*systeme.getObjet(j).getMasse()/norme(distanceVector);
			}
		}
		if(k==0)energie0 = energie;
		energie = (energie-energie0)/energie0*100.0; // Energie devient erreur relative
		energieMeca.push_back(energie);
	}




	// ECRITURE FICHIER
	ofstream f;	
	if(sorties[0]){
	
		// Ecriture des coordonnees dans des fichiers separes
			for(int i = 0 ; i < (int) coordinates.size() ; i++){ // Pour toutes les planètes
			f.open("../positions/"+methode+"_"+systeme.getObjet(i).getLien(),fstream::app); // Récupère le lien
			for(int j = 0 ; j < (int) coordinates[i].size(); j++){ // tous les points
				f << setprecision(20) << (j*h)/(3600*24*365.25) << "\t";
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
		for(int i = 1 ; i < (int) aires.size() ; i++){ // Pour toutes les planètes SAUF SOLEIL
			f.open("../aires/"+methode+"_"+systeme.getObjet(i).getLien(),fstream::app); // Récupère le lien
			for(int j = 0 ; j < (int) aires[i].size(); j++){
				f << setprecision(20) << (j*h)/(3600*24*365.25) << "\t" << aires[i][j] << endl;
			}
			f.close();
		}
	}


	if(sorties[2]){
		// Ecriture des énergie
		f.open("../energies/energieMecanique"+methode +".txt",fstream::app);
		for(int i = 0 ; i < (int) energieMeca.size() ; i++){
			f << setprecision(20) << (i*h)/(3600*24*365.25) << "\t" << energieMeca[i] << endl;
		}
		f.close();
	}

	return systeme;
}


