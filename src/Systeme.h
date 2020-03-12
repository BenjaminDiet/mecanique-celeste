#include <iostream>
#include <vector>
#include <string>

using namespace std;

class Corps;

class Systeme
{
public:
  // Constructeurs
	Systeme(string flux);

  //Accesseurs et mutateurs
	Corps getObjet(int i) const; // Renvoie le Corps i de l'attribut objets
	vector<Corps> getObjets() const;	// Renvoie l'ensemble des objets
	vector<vector <double>> getPositions() const;	// Renvoie l'ensemble des positions des objets
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
	void calculerAires(double h);
private: 
	vector<double> posBarycentre; // Barycentre ?
		
	vector<Corps> objets;
};
