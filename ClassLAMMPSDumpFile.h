#include<fstream>
#include<iostream>
#include<cstdio>
#include<string>
#include<cstdlib>
#include<iomanip>
#include<limits>
#include<math.h>
#include<sstream>


/* Finds the round-off of a number, can be used to find minimum image */
#define ANINT(number) floor(number+0.5)
/*Periodic Boundary Distance*/
#define PBCDistance(Distance,box) Distance-(box)*ANINT(Distance/box)
/*3D Distance Square*/
#define DistanceSquare(a,b,c) (a*a)+(b*b)+(c*c)

using namespace std;
class Coord
{
public:
	double x,y,z;
};
inline double CalculateDistance(Coord &A,Coord &B)
{
	double distance = (A.x - B.x)*(A.x - B.x) + (A.y - B.y)*(A.y - B.y) + (A.z - B.z)*(A.z - B.z);
	return sqrt(distance);
}

inline double MagnitudeVectorXY(Coord &A)
{
	double Mag = (A.x*A.x) + (A.y*A.y);
	return sqrt(Mag);
}

inline double MagnitudeVectorXZ(Coord &A)
{
	double Mag = (A.x*A.x) + (A.z*A.z);
	return sqrt(Mag);
}

class AtomsData
{
public:	
	int AtomID,MoleculeID,AtomTypeID;
	float Charge,mass;
	Coord AtomCoord;
	string temp;
		
	void rd(ifstream &in)
	{
		in>>AtomID>>AtomCoord.x>>AtomCoord.y>>AtomCoord.z>>temp>>temp>>temp>>temp>>temp>>temp;
	}
	void wr(ofstream &out)
	{
		out << AtomID<<" "<<MoleculeID<<" "<<AtomTypeID<<" "<<Charge<<" "<<mass<<" "<<AtomCoord.x<<" "<<AtomCoord.y<<" "<<AtomCoord.z<<endl;
	}
};

// Sphere atom data is based on id, radius xu yu zu (For granular simulations)
class AtomsDataSphere
{
public:	
	int AtomID;
	float Charge,Radius;
	Coord AtomCoord;
	string temp;
		
	void rd(ifstream &in)
	{
		in>>AtomID>>Radius>>AtomCoord.x>>AtomCoord.y>>AtomCoord.z;
	}
	void wr(ofstream &out)
	{
		out << AtomID<<" "<<Radius<<" "<<AtomCoord.x<<" "<<AtomCoord.y<<" "<<AtomCoord.z<<endl;
	}
};

// This class should be used for granular simulations with shear deformation
class DumpFileDataSphere
{
public:
	int TimeStep, NumberofAtoms;
	Coord BoxMin,BoxMax, TiltFac;
	AtomsDataSphere *Atoms;

	void ReadfromFile(string Filename)
	{
		ifstream in(Filename.c_str());
		if(!in)
		{
			cout<<"Error reading "<<Filename;
			cin>>Filename;
		}

		string tempstring;
		in>>tempstring>>tempstring;
		in>>TimeStep;
		in>>tempstring>>tempstring>>tempstring>>tempstring;
		in>>NumberofAtoms;
		in>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring;
		in>>BoxMin.x>>BoxMax.x>>TiltFac.x>>BoxMin.y>>BoxMax.y>>TiltFac.y>>BoxMin.z>>BoxMax.z>>TiltFac.z;
		in>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring;
		Atoms = new AtomsDataSphere[NumberofAtoms];
		for(int AtomsCounter=0;AtomsCounter<NumberofAtoms;++AtomsCounter)
		{
			Atoms[AtomsCounter].rd(in);
		}

		in.close();
	}
	void WritetoFile(string OutFile)
	{
		ofstream out(OutFile.c_str());
		out << "ITEM: TIMESTEP"<<endl;out<<TimeStep<<endl;
		out <<"ITEM: NUMBER OF ATOMS"<<endl;out << NumberofAtoms<<endl;
		out <<"ITEM: BOX BOUNDS xy xz yz pp pp pp"<<endl;out<<BoxMin.x<<" "<<BoxMax.x<<" "<<TiltFac.x<<endl;
												out <<BoxMin.y<<" "<<BoxMax.y<<" "<<TiltFac.y<<endl;
												out <<BoxMin.z<<" "<<BoxMax.z<<" "<<TiltFac.z<<endl;
		out <<"ITEM: ATOMS id radius xu yu zu "<<endl;
		for (int i=0;i < NumberofAtoms;++i)
		{
			Atoms[i].wr(out);
			
		}
		out.close();
	}
			
	void DestroyDumpFileData()
	{
		delete [] Atoms; Atoms=NULL;
	}
};




class DumpFileData
{
public:
	int TimeStep, NumberofAtoms;
	Coord BoxMin,BoxMax;
	AtomsData *Atoms;

	void ReadfromFile(string Filename)
	{
		ifstream in(Filename.c_str());
		if(!in)
		{
			cout<<"Error reading "<<Filename;
			cin>>Filename;
		}

		string tempstring;
		in>>tempstring>>tempstring;
		in>>TimeStep;
		in>>tempstring>>tempstring>>tempstring>>tempstring;
		in>>NumberofAtoms;
		in>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring;
		in>>BoxMin.x>>BoxMax.x>>tempstring>>BoxMin.y>>BoxMax.y>>tempstring>>BoxMin.z>>BoxMax.z>>tempstring;
		in>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring>>tempstring;
		Atoms = new AtomsData[NumberofAtoms];
		for(int AtomsCounter=0;AtomsCounter<NumberofAtoms;++AtomsCounter)
		{
			Atoms[AtomsCounter].rd(in);
		}

		in.close();
	}
	void WritetoFile(string OutFile)
	{
		ofstream out(OutFile.c_str());
		out << "ITEM: TIMESTEP"<<endl;out<<TimeStep<<endl;
		out <<"ITEM: NUMBER OF ATOMS"<<endl;out << NumberofAtoms<<endl;
		out <<"ITEM: BOX BOUNDS pp pp pp"<<endl;out<<BoxMin.x<<" "<<BoxMax.x<<endl;
												out <<BoxMin.y<<" "<<BoxMax.y<<endl;
												out <<BoxMin.z<<" "<<BoxMax.z<<endl;
		out <<"ITEM: ATOMS id mol type q mass xu yu zu "<<endl;
		for (int i=0;i < NumberofAtoms;++i)
		{
			Atoms[i].wr(out);
			
		}
		out.close();
	}
			
	void DestroyDumpFileData()
	{
		delete [] Atoms; Atoms=NULL;
	}
};
string FilenameDotTimestep(string InitialFilename, int TimeStep)
{
	/* Convert Timestep to append to filename */
	stringstream TimeStream;
	TimeStream<<TimeStep;
	return(InitialFilename+'.'+TimeStream.str());
}
void InitializeDumpFile(DumpFileData &File1, DumpFileData &File2)
{
	File2.BoxMax.x = File1.BoxMax.x;
	File2.BoxMax.y = File1.BoxMax.y;
	File2.BoxMax.z = File1.BoxMax.z;
	File2.BoxMin.x = File1.BoxMin.x;
	File2.BoxMin.y = File1.BoxMin.y;
	File2.BoxMin.z = File1.BoxMin.z;
	File2.TimeStep = File1.TimeStep;
}

double Calculate_b(double &a,double &b,double &c)
{
	return a -0.5*(b+c);
}
double Calculate_c(double &a,double &b)
{
	return (a -b);
}
double Calculate_K2 (double &b,double &c, double &Rg2)
{
	return (b*b + (3/4)*c*c)/(Rg2*Rg2);
}
inline double CalculateSquaredDisplacement(Coord A,Coord B)
{
	Coord SquaredDisplacement;
	SquaredDisplacement.x=(A.x-B.x)*(A.x-B.x);
	SquaredDisplacement.y=(A.y-B.y)*(A.y-B.y);
	SquaredDisplacement.z=(A.z-B.z)*(A.z-B.z);

	return((SquaredDisplacement.x+SquaredDisplacement.y+SquaredDisplacement.z));
}

inline double CalculateDisplacement(Coord &A,Coord &B) 
{
	Coord SquaredDisplacement;
	SquaredDisplacement.x= (A.x-B.x)*(A.x-B.x);
	SquaredDisplacement.y=(A.y-B.y)*(A.y-B.y);
	SquaredDisplacement.z= (A.z-B.z)*(A.z-B.z);

	return(sqrt(SquaredDisplacement.y));
}
int FindMin(int a, int b)
{
	int Min = a;
	if (Min > b) Min = b;
	return Min;
}