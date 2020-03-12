#include <iostream>
#include <vector>
#include <string>

using namespace std;

class Corps
{
public:
  // Constructeurs
	Corps();
	Corps(vector<double> r, vector<double> v, vector<double>a, double m, string s);

  //Accesseurs et mutateurs
	vector<double> getPosition() const;
	vector<double> getVitesse() const;
	vector<double> getAcc() const;
	double getMasse() const;
	double getAire() const;
	string getLien() const;

  // Autres méthodes
	void addAcc(vector <double> acce); // Ajout du vecteur acce à l'accélération
	void setAcc(vector <double> acce); // Ajout du vecteur acce à l'accélération
	void saveAcc();	// Enregistre l'accélération dans accAvant
	void emptyAcc();	// Mets l'accélération au vecteur nul
	void loiDesAires(vector <double> barycentre, double h);

	
  // Résolution
	void majPositionVerlet(double h);
	void majVitesseVerlet(double h);

	void majPositionEuler(double h);
	void majVitesseEuler(double h);


	// RK4	
	void AddPosition(vector<double> k,double h,double indice);
	void SetPosition(vector<double> k,double indice);
	void SubPosition(vector<double> k,double h,double indice);
	void SetVitesse(vector<double> k,double indice);


private: 
	vector<double> position;
	vector<double> vitesse;
	vector<double> acc;

	vector<double> accAvant; // Stockage de l'accélération pour Verlet

	double masse;

	double aire;

	string lien;
};
