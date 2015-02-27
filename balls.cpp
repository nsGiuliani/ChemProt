#include<iostream>
#include<cmath>
#include<math.h>
#include<cstdlib>//needed for rand
#include<ctime> //needed for srand(time(0))
#include<fstream> //needed for files
#include<vector>  
#include "atom.cpp"

using namespace std;

int natoms = 5;

int main(){
	ofstream fout;
	  
	fout.open("5atom.xyz");

	fout << natoms << endl;
	fout << "comment line" << endl;

	float x1 = 0.5*50;
	float y1 = 0.25*50;
	float z1 = 1*50;
	float x2 = 0.5*50;
	float y2 = 0.5*50;
	float z2 = 1*50;
	float x3 = 0.75*50;
	float y3 = 0.5*50;
	float z3 = 1*50;
	float x4 = 0.75*50;
	float y4 = 0.75*50;
	float z4 = 1*50;
	float x5 = 0.5*50;
	float y5 = 0.75*50;
	float z5 = 1*50;
	vector<Atom> atoms;

	int iterations = 10000;
	double box = 100;
	double timeStep = .0000000000000001;
	double ang = 0.0000000001;
	double x;
	double y;
	double z;
	int numAtoms = 5;
	double E = 1.712*(.000000000000000000001);
	double sig = 3.4*(.0000000001);
	double Fx;
	double Fy;
	double Fz;
	double dist;
	double dx=0.0;
	double dy=0.0;
	double dz=0.0;

	//temp constants for Hooke's law
	double hX = 10;
	double hK = .000001;
	double hF;

	atoms.push_back(Atom(39.948*1.66053892 * pow(10,-27)*3, x1, y1, z1, 0, 0, 0));
	fout << "Ar" << " ";
	fout << x1 << " " << y1 << " " << z1 << endl;
		
	atoms.push_back(Atom(39.948*1.66053892 * pow(10,-27)*3, x2, y2, z2, 0, 0, 0));
	fout << "Ar" << " ";
	fout << x2 << " " << y2 << " " << z2 << endl;
		
	atoms.push_back(Atom(39.948*1.66053892 * pow(10,-27)*3, x3, y3, z3, 0, 0, 0));
	fout << "Ar" << " ";
	fout << x3 << " " << y3 << " " << z3 << endl;

	atoms.push_back(Atom(39.948*1.66053892 * pow(10,-27)*3, x4, y4, z4, 0, 0, 0));
	fout << "Ar" << " ";
	fout << x4 << " " << y4 << " " << z4 << endl;

	atoms.push_back(Atom(39.948*1.66053892 * pow(10,-27)*3, x5, y5, z5, 0, 0, 0));
	fout << "Ar" << " ";
	fout << x5 << " " << y5 << " " << z5 << endl;

	for (int i=0; i<iterations; i++ ) {
		fout << numAtoms << endl;
		fout << "iteration " << i << endl;
		for (int k=0; k<numAtoms; k++) {
			fout<< "Ar ";
			fout<< atoms[k].getxCoor() <<"	"<<atoms[k].getyCoor() << "	"<< atoms[k].getzCoor()<< endl;
		}
		for (int k=0; k<numAtoms; k++) { //k is the atom which is being updated
			//update the coordinates and make sure they are in the box
			x= atoms[k].getxCoor() + atoms[k].getxVel()*timeStep/ang;
			y= atoms[k].getyCoor() + atoms[k].getyVel()*timeStep/ang;
			z= atoms[k].getzCoor() + atoms[k].getzVel()*timeStep/ang;
			if(x >= box){
				x =fmod(x, box);
			}
			if(x <0) {
				x =fmod(x, box) + box;
			}
			if(y >= box){
				y = fmod(y, box);
			}
			if(y <0) {
				y =fmod(y, box) +box;
			}
			if(z >= box){
				z =fmod(z, box);
			}
			if(z <0) {
				z =fmod(z, box) + box;
			}

			atoms[k].set_Coor(x,y,z);
		}
		for (int k=0; k<numAtoms; k++) { //k is the atom which is being updated
			for (int l = 0; l<numAtoms;l++){ //l is the atom whose presence is applying force to k
				double min = 1000000000;
				if(l != k){ //make sure they are not the same atom
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());
					if (dist < min) {
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor()-box)*(atoms[l].getxCoor()-atoms[k].getxCoor()-box)
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (dist < min) {
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor()-box;
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}			
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor()+box)*(atoms[l].getxCoor()-atoms[k].getxCoor()+box)
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (dist < min){
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor()+box;
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}			
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor()-box)*(atoms[l].getyCoor()-atoms[k].getyCoor()-box)
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (dist < min){
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor()-box;
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor()+box)*(atoms[l].getyCoor()-atoms[k].getyCoor()+box)
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (dist < min){
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor()+box;
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}		
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor()-box)*(atoms[l].getzCoor()-atoms[k].getzCoor()-box);					
					if (dist < min){
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor()-box;
					}
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor()+box)*(atoms[l].getzCoor()-atoms[k].getzCoor()+box);					
					if (dist < min){
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor()+box;
					}
					//Hooke's portion
					min = sqrt(min);
					if (k == 0){
						if(l ==1){
							hF = (-1)*(hX- abs(min))*hK;
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
							hF = (-1)*(hX- abs(min))*hK;
							Fx = hF * (dx/min);
							Fy = hF * (dy/min);
							Fz = hF * (dz/min);
							x= atoms[k].getxVel() + (Fx/atoms[k].getMass())*timeStep;
							y= atoms[k].getyVel() + (Fy/atoms[k].getMass())*timeStep;
							z= atoms[k].getzVel() + (Fz/atoms[k].getMass())*timeStep;
							atoms[k].set_Vel(x,y,z);
							cout <<min << " " << dx<< " " << Fx << endl;
						}
					}
					if (k != 0 && k != numAtoms -1){
						if(l == k-1 || l == k+1){
							hF = (-1)*(hX- abs(min))*hK;
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
    



