#include "Systeme.h"
#include "Corps.h"
#include "Utiles.h"
#include "Constantes.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;


Systeme::Systeme(string flux,double h):posBarycentre(3){
	
	// LECTURE DES DONNEES
	fstream input("../data/" + flux, ios_base::in);
	vector <vector <double>> positions; 
	vector <vector <double>> vitesses;
	vector <vector <double>> accelerations;
	vector <double> masses;	
	vector<string> liens;	
	vector<string> noms;
	vector<int> natures;		
	vector<double> tmpPlusH(3);
	vector<double> tmpMoinsH(3);
	vector<double> VitBarycentre(3);

	double rX, rY, rZ, vX, vY, vZ, m;	// Valeurs temporaires
	string l, n;					// Idem
	int nat;

	while(!input.eof()){
		vector<double> r, v, a;
		input >> n >> nat >> rX >> rY >> rZ >> vX >> vY >> vZ >> m;
		r.push_back(rX); r.push_back(rY); r.push_back(rZ);
		v.push_back(vX); v.push_back(vY); v.push_back(vZ);
		a.push_back(0); a.push_back(0); a.push_back(0);
		positions.push_back(r);
		vitesses.push_back(v);
		accelerations.push_back(a);
		masses.push_back(m);
		natures.push_back(nat);
		
		
		liens.push_back("data"+n+".txt");	
		noms.push_back(n);		
	}
	input.close();
	// Nettoie les fins de fichiers
	positions.pop_back();
	vitesses.pop_back();
	accelerations.pop_back();
	masses.pop_back();
	liens.pop_back();
	// ALLOUCHE : supprimer ici la vitesse du centre de masse du systeme. On calcule la vitesse du centre de masse puis on soustrait
	// cette vitesse à chaque astre 

	// INITIALISATION DES CORPS
	for(size_t i = 0 ; i < positions.size() ; i++){
		objets.push_back(Corps(positions[i], vitesses[i], accelerations[i], masses[i], liens[i], noms[i], natures[i]));
	}


	// BARYCENTRE
	this->calculerBarycentre();
	this->centrerBarycentre();

	this->resoudreVerlet(-h, false);
	this->calculerBarycentre();
	tmpMoinsH = posBarycentre;

	this->resoudreVerlet(h, false);
	this->resoudreVerlet(h, false);

	this->calculerBarycentre();
	tmpPlusH = posBarycentre;

	this->resoudreVerlet(-h, false);

	for(int i = 0 ; i < 3 ; i++){
		VitBarycentre[i] = (1.0/(2.0*h))*(tmpPlusH[i] - tmpMoinsH[i]);
	}
	
	for(size_t i = 0 ; i < objets.size() ; i++){
		objets[i].SetVitesse(objets[i].getVitesse()-VitBarycentre, 0); // Soustrait vitesse du barycentre
	}
	
}


Systeme &Systeme::operator=(const Systeme &source){
	this->posBarycentre=source.posBarycentre;
	this->objets=source.getObjets();
			
	return *this;
}




// ACCESSEURS
vector<Corps> Systeme::getObjets() const{return this->objets;}

Corps &Systeme::operator[](size_t index){
	if(index >= objets.size()){
		cerr << "Index hors limites" << endl;
		return objets[0];
	}
	return objets[index];
}



vector<vector <double>> Systeme::getPositions() const{
	vector<vector <double>> pos;
	for(size_t i = 0 ; i < objets.size() ; i++){
		pos.push_back(objets[i].getPosition());
	}
	return pos;
}



vector<vector <double>> Systeme::getVitesses() const{
	vector<vector <double>> vit;
	for(size_t i = 0 ; i < objets.size() ; i++){
		vit.push_back(objets[i].getVitesse());
	}
	return vit;
}


vector<double> Systeme::getAires() const{
	vector<double>aires;
	for(size_t i = 0 ; i < this->objets.size() ; i++){
		aires.push_back(objets[i].getAire());
	}
	return aires;
}


int Systeme::getSize() const{return this->objets.size();}



// Calculs
void Systeme::calculerBarycentre()
{
		double massetot(0);
		
		for (size_t i=0;i<objets.size();i++){
				for(int j = 0 ; j < 3 ; j++){posBarycentre[j]+=objets[i].getMasse()*objets[i].getPosition()[j];}
				massetot+=objets[i].getMasse();
		}
		
		for(int i = 0 ; i < 3 ; i++){posBarycentre[i]/=massetot;}		
}

void Systeme::centrerBarycentre(){	
		for(size_t i=0; i < objets.size() ; i++){objets[i].SubPosition(posBarycentre,1,0);} // Enlève coordonnées barycentre
}






// RESOLUTION

void Systeme::resoudreEuler(double h, bool relativiste){


	// On calcule a(t+dt)
	if(relativiste) this->calculerAccRelativite();
	else this->calculerAcc();
	for(size_t i = 0 ; i < this->objets.size() ; i++){	
		this->objets[i].majPositionEuler(h); // On en déduit la position
		this->objets[i].majVitesseEuler(h); // On en déduit la vitesse
	}

	
}

void Systeme::resoudreEulerCromer(double h, bool relativiste){

	// On calcule a(t+dt)
	if(relativiste) this->calculerAccRelativite();
	else this->calculerAcc();
	for(size_t i = 0 ; i < this->objets.size() ; i++){	
		this->objets[i].majVitesseEuler(h); // On en déduit la position
		this->objets[i].majPositionEuler(h); // On en déduit la position
	}



}


void Systeme::resoudreVerlet(double h, bool relativiste){

	for(size_t i = 0 ; i < this->objets.size() ; i++){		
		this->objets[i].majPositionVerlet(h); // On calcule x(t+dt)
	}
	
	for(size_t i = 0 ; i < this->objets.size() ; i++){		
		this->objets[i].majVitesseVerlet(h); // On calcule v(t+dt)
	}
	if(relativiste) this->calculerAccRelativite();
	else this->calculerAcc();

	for(size_t i = 0 ; i < this->objets.size() ; i++){		
		this->objets[i].majVitesseVerlet(h); // On calcule v(t+dt)
	}
}


	

void Systeme::calculerAcc(){

	for(int i = 0 ; i < (int) this->objets.size() ; i++){
		this->objets[i].emptyAcc(); // On vide l'accélération
	}
	for(int i = 0 ; i < (int) this->objets.size() ; i++){
	

		for(int j = i+1 ; j < (int) this->objets.size() ; j++){
			
			vector<double> vecteurIJ(3),  accI(3), accJ(3);	
	
			vecteurIJ = distance(this->objets[i].getPosition(), this->objets[j].getPosition());
			

			
			double normeVecteur = norme(vecteurIJ);
			double normeCarree = normeVecteur*normeVecteur;
			
			double forceI = this->objets[j].getMasse()*G/normeCarree;
			double forceJ = this->objets[i].getMasse()*G/normeCarree;


			for(int k = 0 ; k < (int) vecteurIJ.size() ; k++){
				accI[k] = vecteurIJ[k]*forceI/normeVecteur;
				accJ[k] = -vecteurIJ[k]*forceJ/normeVecteur;				
			}
			
			this->objets[i].addAcc(accI);
			this->objets[j].addAcc(accJ);
			 
		} 
	}

}



void Systeme::calculerAccRelativite(){
	


	calculerAcc();
	

	for(int i = 0 ; i < (int) this->objets.size() ; i++){

		vector <double> a = this->objets[i].getAcc();
		double alpha(0);
		double beta(0);
		double delta(0);
		double gamma(0);

		for(int j = 0 ; j < (int) this->objets.size() ; j++){if(i!=j){
			vector <double> ab = this->objets[j].getAcc();
			vector <double> va = this->objets[i].getVitesse();
			vector <double> vb = this->objets[j].getVitesse();
			vector <double> rba = distance(this->objets[i].getPosition(), this->objets[j].getPosition());
			vector <double> rab = oppose(rba);
			vector <double> nba = normaliser(rba);
			vector <double> nab = normaliser(rab);
			double mb = this->objets[j].getMasse();

			delta = G*mb/pow(norme(rab), 2)/pow(c,2);
			beta = delta * produitScalaire(nab, distance(multScalaire(3, vb),multScalaire(4, va)));

			double epsilon(0);
			for(int k = 0 ; k < (int) this->objets.size() ; k++){
					double mc = this->objets[k].getMasse();
					vector <double> rac = distance(this->objets[k].getPosition(), this->objets[i].getPosition());
					vector <double> rbc = distance(this->objets[k].getPosition(), this->objets[j].getPosition());
				if(i!=k) epsilon -= 4.0*G*mc/norme(rac);
				if(j!=k) epsilon -= G*mc/norme(rbc);
			}
	
			alpha = delta*(epsilon + pow(norme(va),2) +2.0*pow(norme(vb),2) - 4.0*produitScalaire(va, vb) -1.5*pow(produitScalaire(nab,vb),2)+0.5*produitScalaire(rba,ab));
			gamma = G*mb/norme(rab)*3.5/pow(c,2);


			a = add(a, multScalaire(alpha, nba));
			a = add(a, multScalaire(beta, distance(vb, va)));
			a = add(a, multScalaire(gamma,ab));
		}}		
		
		this->objets[i].setAcc(a);
	}
	

}




// RK4
void Systeme::resoudreRK4(double h, bool relativiste){
	
	
		vector<double> kr1(objets.size()*3);
		vector<double> kr2(objets.size()*3);
		vector<double> kr3(objets.size()*3);
		vector<double> kr4(objets.size()*3);
		
		vector<double> kv1(objets.size()*3);
		vector<double> kv2(objets.size()*3);
		vector<double> kv3(objets.size()*3);
		vector<double> kv4(objets.size()*3);
		int compteur(0);
		int indice(0);
		
		if(relativiste) this->calculerAccRelativite();
	else this->calculerAcc();
		
		for (int i=0;i<(int) objets.size();i++)
		{	 
			for(int j=0;j<3;j++)//calcul du coeff k1 
			{		
					kr1[compteur]=objets[i].getVitesse()[j]*h;//pour la position k1=V1*h
					kv1[compteur]=objets[i].getAcc()[j]*h;//pour la vitesse K1=a(r)
					compteur++;//permet de mettre les k1 de chaque corps dans un seul vector
			}
			objets[i].emptyAcc();//dans la méthode CalculerAcc() on ajoute les accélérations au vecteur accélération donc il faut remettre ce vecteur à zero quand on veut calculer une nouvelle accélération 
		}
		compteur=0;
		
		for(int i=0;i<(int)objets.size();i++)//on met a jour la position pour calculer a(ri+k1*0.5)
		{		
				objets[i].AddPosition(kr1,0.5,indice);
				indice+=3;// permet de savoir quand commence les indice du prochain objet dans le vecteur k qui contient tous les indices de tous les objets à la suite.
		}
		indice=0;
		
		if(relativiste) this->calculerAccRelativite();
	else this->calculerAcc();//on calcule a(ri+k1*0.5)
		
		for(int i=0;i<(int)objets.size();i++)//on remet la position comme avant le calcul de l'accélération
		{	
				objets[i].SubPosition(kr1,0.5,indice);
				indice+=3;
		}
	
		indice=0;
		
		for (int i=0;i<(int) objets.size();i++)//calcul du coeff k2
		{
			for(int j=0;j<3;j++)
			{
					kr2[compteur]=(objets[i].getVitesse()[j]+0.5*kv1[compteur])*h;//pour la position k2=(V+0.5*kV1)*h
					kv2[compteur]=objets[i].getAcc()[j]*h;//pour la vitesse k2=a(r+kr1*0.5)*h
					compteur++;	
			}
			objets[i].emptyAcc();
		}
		compteur=0;
		
		for(int j=0;j<(int)objets.size();j++)//on met a jour la position pour calculer a(ri+k2*0.5)
		{		
				objets[j].AddPosition(kr2,0.5,indice);
				indice+=3;
		}
		indice=0;
		if(relativiste) this->calculerAccRelativite();
	else this->calculerAcc();//on calcule a(ri+k2*0.5)
		
		for(int j=0;j<(int)objets.size();j++)//on remet la position comme avant le calcul de l'accélration
		{
				objets[j].SubPosition(kr2,0.5,indice);
				indice+=3;
		}
		indice=0;
		
		for (int i=0;i<(int) objets.size();i++)//calcul coeff k3
		{
			for(int j=0;j<3;j++)
			{
					kr3[compteur]=(objets[i].getVitesse()[j]+0.5*kv2[compteur])*h;//pour la position k3=(V+0.5*kV2)*h
					kv3[compteur]=objets[i].getAcc()[j]*h;//pour la vitesse k3=a(r+kr2*0.5)*h
					compteur++;	
			}
			objets[i].emptyAcc();
		}
		compteur=0;
		
		for(int j=0;j<(int)objets.size();j++)//on met a jour la position pour calculer a(ri+k3)
		{		
			
				objets[j].AddPosition(kr3,1,indice);
				indice+=3;
			
		}
		indice=0;
		
		if(relativiste) this->calculerAccRelativite();
	else this->calculerAcc();//on calcule a(ri+k3)
		
		for(int j=0;j<(int)objets.size();j++)//on remet la position comme avant le calcul de l'accélération
		{
				objets[j].SubPosition(kr3,1,indice);
				indice+=3;
		}
		indice=0;
		
		for (int i=0;i<(int) objets.size();i++)//calcul coeff k4
		{
			for(int j=0;j<3;j++)
			{
					kr4[compteur]=(objets[i].getVitesse()[j]+kv3[compteur])*h;//pour la position k4=(V+kV3)*h
					kv4[compteur]=objets[i].getAcc()[j]*h;//pour la vitesse k4=a(r+kr3)*h
					compteur++;
			}
			objets[i].emptyAcc();
		}
		compteur=0;
		vector<double> pos(kr1.size());
		vector<double> vit(kv1.size());
		for (int i=0;i<(int) objets.size();i++)
		{
			for(int j=0;j<3;j++)
				{
					pos[compteur]=objets[i].getPosition()[j]+1.0/6.0*(kr1[compteur]+2*kr2[compteur]+2*kr3[compteur]+kr4[compteur]);// on calcule la position à l'indice n+1
					vit[compteur]=objets[i].getVitesse()[j]+1.0/6.0*(kv1[compteur]+2*kv2[compteur]+2*kv3[compteur]+kv4[compteur]);//on calcule la vitesse à l'indice n+1
					compteur++;
				}	
					objets[i].SetPosition(pos,indice);// on remplace l'ancienne position par la nouvelle 
					objets[i].SetVitesse(vit,indice);//on remplace l'ancienne vitesse par la nouvelle
					indice+=3;
					
		}	
}
