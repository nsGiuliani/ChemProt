#include<iostream>
#include<cmath>
#include<math.h>
#include<cstdlib>//needed for rand
#include<ctime> //needed for srand(time(0))
#include<fstream> //needed for files
#include<vector>  
#include "acid.cpp"
#include <string>

using namespace std;


int main(){

	vector<Acid> atoms;

	int iterations = 10000;
	double timeStep = .0000000000000001;
	double ang = 0.0000000001;
	double x;
	double y;
	double z;
	double E = 1.712*(.000000000000000000001);
	double sig = 3.4*(.0000000001);
	double Fx;
	double Fy;
	double Fz;
	double dist;
	double dx=0.0;
	double dy=0.0;
	double dz=0.0;
	int numAtoms;
	string comment;

	double c1,c2,c3;
	string name;
	char mass[10];
	

	//temp constants for Hooke's law
	vector<double> hX;
	double hK = 5000 * ang;
	double hF;

	ifstream prot_file;
	prot_file.open("dumb.txt", ios::in);
	prot_file >> numAtoms;

	prot_file >> comment;

	for(int i = 0; i < numAtoms; i++) {
		prot_file >> name;
		prot_file >> mass;
		prot_file >> c1;
		prot_file >> c2;
		prot_file >> c3;
		atoms.push_back(Acid(name, atof(mass), c1,c2,c3, 0, 0, 0));
		cout <<name<<mass<< c1<<c2<<c3<< endl;

	}
	prot_file.close();	

	//set the desired distances and save them in a vector

	for (int i=0; i<numAtoms-1; i++) {
		hX.push_back(sqrt((atoms[i].getxCoor()-atoms[i+1].getxCoor())*(atoms[i].getxCoor()-atoms[i+1].getxCoor())
				+(atoms[i].getyCoor()-atoms[i+1].getyCoor())*(atoms[i].getyCoor()-atoms[i+1].getyCoor())
				+ (atoms[i].getzCoor()-atoms[i+1].getzCoor())*(atoms[i].getzCoor()-atoms[i+1].getzCoor())));
	}

	ofstream fout;
	  
	fout.open("trp.xyz");

	for (int i=0; i<iterations; i++ ) {
		fout << numAtoms << endl;
		fout << "iteration " << i << endl;
		for (int k=0; k<numAtoms; k++) {
			fout<<atoms[k].getName()<< "  " << atoms[k].getxCoor() <<"	"<<atoms[k].getyCoor() << "	"<< atoms[k].getzCoor()<< endl;
		}
		for (int k=0; k<numAtoms; k++) { //k is the atom which is being updated
			//update the coordinates and make sure they are in the box
			x= atoms[k].getxCoor() + atoms[k].getxVel()*timeStep/ang;
			y= atoms[k].getyCoor() + atoms[k].getyVel()*timeStep/ang;
			z= atoms[k].getzCoor() + atoms[k].getzVel()*timeStep/ang;
			

			atoms[k].set_Coor(x,y,z);
		}
		for (int k=0; k<numAtoms; k++) { //k is the atom which is being updated
			for (int l = 0; l<numAtoms;l++){ //l is the atom whose presence is applying force to k
				double min = 1000000000;
				if(l != k){ //make sure they are not the same atom
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());
					
					min = dist;
					dx =atoms[l].getxCoor()-atoms[k].getxCoor();
					dy =atoms[l].getyCoor()-atoms[k].getyCoor();
					dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					
					
					//Hooke's portion
					min = sqrt(min);
					hK = 5.552 * pow(10,28) * ((atoms[k].getMass()*atoms[l].getMass()) / (atoms[k].getMass()+atoms[k].getMass())) * ang;
					if (k == 0){
						if(l ==1){
							hF = (-1)*(hX[k]- abs(min))*hK;
							Fx = hF * (dx/min);
							Fy = hF * (dy/min);
							Fz = hF * (dz/min);
							x= atoms[k].getxVel() + (Fx/atoms[k].getMass())*timeStep;
							y= atoms[k].getyVel() + (Fy/atoms[k].getMass())*timeStep;
							z= atoms[k].getzVel() + (Fz/atoms[k].getMass())*timeStep;
							atoms[k].set_Vel(x,y,z);
						}
					}
					if (k == numAtoms -1){
						if(l == k-1){
							hF = (-1)*(hX[k-1]- abs(min))*hK;
							Fx = hF * (dx/min);
							Fy = hF * (dy/min);
							Fz = hF * (dz/min);
							x= atoms[k].getxVel() + (Fx/atoms[k].getMass())*timeStep;
							y= atoms[k].getyVel() + (Fy/atoms[k].getMass())*timeStep;
							z= atoms[k].getzVel() + (Fz/atoms[k].getMass())*timeStep;
							atoms[k].set_Vel(x,y,z);
						}
					}
					if (k != 0 && k != numAtoms -1){
						if(l == k-1){
							hF = (-1)*(hX[k-1]- abs(min))*hK;
							Fx = hF * (dx/min);
							Fy = hF * (dy/min);
							Fz = hF * (dz/min);
							x= atoms[k].getxVel() + (Fx/atoms[k].getMass())*timeStep;
							y= atoms[k].getyVel() + (Fy/atoms[k].getMass())*timeStep;
							z= atoms[k].getzVel() + (Fz/atoms[k].getMass())*timeStep;
							atoms[k].set_Vel(x,y,z);
						}
						if(l == k+1){
							hF = (-1)*(hX[k]- abs(min))*hK;
							Fx = hF * (dx/min);
							Fy = hF * (dy/min);
							Fz = hF * (dz/min);
							x= atoms[k].getxVel() + (Fx/atoms[k].getMass())*timeStep;
							y= atoms[k].getyVel() + (Fy/atoms[k].getMass())*timeStep;
							z= atoms[k].getzVel() + (Fz/atoms[k].getMass())*timeStep;
							atoms[k].set_Vel(x,y,z);
						}
					}

					dx= dx*ang;
					dy=dy*ang;
					dz=dz*ang;
					// update the velocity 
					double D = dx*dx+dy*dy+dz*dz;
					double D_7= (dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)
						*(dx*dx+dy*dy+dz*dz);
					double D_4 = (dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz);
					double sig_6= sig*sig*sig*sig*sig*sig;
					Fx = -1*(24*E*sig_6*dx*((2*sig_6/(D_7))-(1/(D_4))));
					Fy = -1*(24*E*sig_6*dy*((2*sig_6/(D_7))-(1/(D_4))));					
					Fz = -1*(24*E*sig_6*dz*((2*sig_6/(D_7))-(1/(D_4))));
					x= atoms[k].getxVel() + (Fx/atoms[k].getMass())*timeStep;
					y= atoms[k].getyVel() + (Fy/atoms[k].getMass())*timeStep;
					z= atoms[k].getzVel() + (Fz/atoms[k].getMass())*timeStep;
					atoms[k].set_Vel(x,y,z);
					//cout << Fx << endl;

				}		
																	
			}
		}
	}

	fout.close();

	

	return 0;

}  
    



