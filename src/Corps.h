#include <iostream>
#include <vector>
#include <string>
#include "parametre.h"

using namespace std;

class Corps :public parametre
{
public:
  // Constructeurs
	Corps();
	Corps(vector<double> r, vector<double> v, vector<double>a, double m, string s, string n, int nature);

  //Accesseurs et mutateurs
	vector<double> getPosition() const;
	vector<double> getVitesse() const;
	vector<double> getAcc() const;
	double getMasse() const;
	double getAire() const;
	string getLien() const;
	int getNature() const;
	string getNom() const;
	double getExcentricite() const;
	double getPeriode() const;

  // Autres méthodes
	void addAcc(vector <double> acce); // Ajout du vecteur acce à l'accélération
	void setAcc(vector <double> acce); // Ajout du vecteur acce à l'accélération
	void multiplierVitesse(double k);
	void emptyAcc();	// Mets l'accélération au vecteur nul
	void setExcentricite(double e);
	void setPeriode(double p);

	
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

	double masse;


	string lien;
	string nom;
	int nature;
	
	
	double aire;
	double periode;
	double excentricite;
};

