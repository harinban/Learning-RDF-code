// Written by Fardin Khabaz (FK)

#include "ClassLAMMPSDumpFile.h"
#include "CalculateRDF_Dump.h"


int main()
{
	int StartTime,
		EndTime,
		IntervalTime,
		Flag;
	string TempString, //Temporary string dumping text from a file
		DumpFilename1,
		InstantaneousResults;
	// Configuration File
	double Lx, Ly, Lz, count, LayerThickness, Rho;

	ifstream ReadConfig("Config.txt"); //Reading from a file
	ReadConfig >> TempString >> DumpFilename1 >> TempString >> StartTime >> TempString >> 
		EndTime >> TempString >> IntervalTime >> TempString >> Flag >> TempString >> LayerThickness;
	//Reading it from ReadConfig

	ReadConfig.close();  //Closing the ReadConfig.

	int NumberofIntervals=int(floor((EndTime-StartTime)/IntervalTime)) + 1; //floor: rounded downard

	
	DumpFileDataSphere D1,D2;
	D1.ReadfromFile(FilenameDotTimestep(DumpFilename1,StartTime).c_str());
	Lx = D1.BoxMax.x - D1.BoxMin.x; 
        Ly = D1.BoxMax.y - D1.BoxMin.y;
	Lz = D1.BoxMax.z - D1.BoxMin.z;
	D1.DestroyDumpFileData();
	double L_half = Lz/2.0;
	double BinSize  = 0.1;
	int NumberofBins = int(floor(L_half/BinSize));
	double *rdf; rdf = new double[NumberofBins];

	double *Ave_gr; Ave_gr= new double[NumberofBins];
	
	double **gr;
	gr = new double *[NumberofBins];
	for (int i = 0; i < NumberofBins; i++)
	{
		gr[i] = new double[NumberofIntervals];
	}
	for (int i = 0; i < NumberofBins; i++)
	{
		for (int j = 0; j < NumberofIntervals; j++) gr[i][j] = 0.0;
	}

	string Filename1, Filename2;
	Filename1 = "PairDist_3D.txt";
	Filename2 = "PairDist_3D_Ave.txt";

	ofstream Out(Filename1.c_str());
	ofstream OutAve(Filename2.c_str());

	for (int i = 0;  i < NumberofBins; i++) 
	{
		rdf[i] = 0.0;
		Ave_gr[i] = 0.0;
	}
	Rho = 0;
	for (int FileIndex = 0; FileIndex<NumberofIntervals; FileIndex++)
	{
		D1.ReadfromFile(FilenameDotTimestep(DumpFilename1, StartTime + FileIndex*IntervalTime).c_str());
		
				RDF_3D(D1, rdf, NumberofBins, BinSize, LayerThickness, Rho);
				for (int i = 0; i < NumberofBins; i++)
				{
					Ave_gr[i] += rdf[i];
					gr[i][FileIndex] = rdf[i];
					rdf[i] = 0;
				}
		D1.DestroyDumpFileData();
	}
	cout <<Rho;
	for (int i = 0; i < NumberofBins; i++)
	{
		double r1 = i*BinSize;
		double r2 = r1 + BinSize;
		Ave_gr[i] /= NumberofIntervals;
		OutAve << (r1+r2)/2.0 << " "<< Ave_gr[i] <<endl;
	}

	for (int i = 0; i < NumberofBins; i++)
	{
		double r1 = i*BinSize;
		double r2 = r1 + BinSize;
		Out << (r1 + r2) / 2.0;
		for (int j = 0; j < NumberofIntervals; j++)
		{
			Out<< " "<<gr[i][j];
		}
		Out << endl;
	}

	//delete [] rdf; rdf = NULL;
	return 9111989;
}

