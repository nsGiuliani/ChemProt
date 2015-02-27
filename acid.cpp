//basic atom class

using namespace std;

class Acid {
	string name;
	double xCoor, yCoor, zCoor, xVel, yVel, zVel, mass;

  public:
	Acid ();
	Acid(string, double, double, double, double, double, double, double);
	void set_Acid(string, double, double, double, double, double, double, double);
	void set_Coor (double, double, double);
	void set_Vel (double, double, double);
	string getName() {
		return name;
	}
	double getMass() {
		return mass;
	}
	double getxCoor() {
		return xCoor;
	}
	double getyCoor() {
		return yCoor;
	}
	double getzCoor() {
		return zCoor;
	}
	double getxVel() {
		return xVel;
	}
	double getyVel() {
		return yVel;
	}
	double getzVel() {
		return zVel;
	}
};
Acid::Acid () {
	name = "?";
	xCoor = 0;
	yCoor = 0;
	zCoor = 0;
	xVel = 0;
	yVel = 0;
	zVel = 0;
	mass = 0;
}

Acid::Acid (string n, double m,double x1, double y1, double z1, double x2, double y2, double z2) {
	xCoor = x1;
	yCoor = y1;
	zCoor = z1;
	xVel = x2;
	yVel = y2;
	zVel = z2;
	mass = m;
	name = n;
}

void Acid::set_Acid (string n, double m,double x1, double y1, double z1, double x2, double y2, double z2) {
	xCoor = x1;
	yCoor = y1;
	zCoor = z1;
	xVel = x2;
	yVel = y2;
	zVel = z2;
	mass = m;
	name = n;
}
void Acid::set_Coor (double x1, double y1, double z1) {
	xCoor = x1;
	yCoor = y1;
	zCoor = z1;
}
void Acid::set_Vel (double x1, double y1, double z1) {
	xVel = x1;
	yVel = y1;
	zVel = z1;
}

