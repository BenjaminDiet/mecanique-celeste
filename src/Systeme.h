#include <iostream>
#include <vector>
#include <string>
#include "parametre.h"

using namespace std;

class Corps;

class Systeme: public parametre
{
public:
  // Constructeurs
	Systeme(string flux,double h);

  //Accesseurs et mutateurs
	Corps &operator[](size_t); // Surcharge l'opérateur indiciel
	Systeme &operator=(const Systeme &source); // Opérateur copie
	vector<Corps> getObjets() const;	// Renvoie l'ensemble des 
	vector<vector <double>> getPositions() const;	// Renvoie l'ensemble des positions des objets
	vector<vector <double>> getVitesses() const;	// Renvoie l'ensemble des vitesses des objets
	vector<double> getAires() const;	// Renvoie l'ensemble des positions des objets
	int getSize() const;	// Renvoie le nombre d'objets

  // Autres méthodes
	void resoudreEuler(double h, bool relativiste);
	void resoudreEulerCromer(double h, bool relativiste);
	void resoudreVerlet(double h, bool relativiste);
	void resoudreRK4(double h, bool relativiste);


	void calculerAcc();	// Calcule l'accélération depuis la force gravitationnelle de tous les corps
	void calculerAccRelativite();	// Calcule l'accélération depuis la force gravitationnelle de tous les corps RELAT

	void calculerBarycentre();
	void centrerBarycentre();
private: 
	vector<double> posBarycentre; // Barycentre ?
		
	vector<Corps> objets;
};

