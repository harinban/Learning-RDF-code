void RDF_3D(DumpFileDataSphere D1, double * rdf, int &NumberofBins, double &BinSize, double &LayerThickness, double &Rho)
{
	double Lx, Ly, Lz,  Distance;
	Coord Delta;
	Lx = D1.BoxMax.y - D1.BoxMin.y; 
    	Ly = D1.BoxMax.y - D1.BoxMin.y;
	Lz = D1.BoxMax.z - D1.BoxMin.z;
	double L_half = Ly/2.0;
	int BinIndex;

	for (int i = 0;  i < NumberofBins; i++) rdf[i] = 0.0;

	double Count = 0.0;
	for (int i = 0; i < D1.NumberofAtoms; i++)
	{
		D1.Atoms[i].AtomCoord.x = PBCDistance(D1.Atoms[i].AtomCoord.x,Lx);
		D1.Atoms[i].AtomCoord.y = PBCDistance(D1.Atoms[i].AtomCoord.y,Ly);
		D1.Atoms[i].AtomCoord.z = PBCDistance(D1.Atoms[i].AtomCoord.z,Lz);

		for (int j = i + 1; j < D1.NumberofAtoms; j++)
		{
			Delta.x = D1.Atoms[i].AtomCoord.x - D1.Atoms[j].AtomCoord.x;
			Delta.x = PBCDistance(Delta.x,Lx);
			Delta.y = D1.Atoms[i].AtomCoord.y - D1.Atoms[j].AtomCoord.y;
			if (Delta.y < 0){
				while (Delta.y < 0.0){
					Delta.y += Ly;
					Delta.x += D1.TiltFac.x;
				}
				if (Delta.y > Ly/2.0){
					Delta.y -= Ly;
					Delta.x -= D1.TiltFac.x;
				}
			}else{
				while (Delta.y > 0.0){
					Delta.y -= Ly;
					Delta.x -= D1.TiltFac.x;
				}
				if (Delta.y < -Ly/2.0){
					Delta.y += Ly;
					Delta.x += D1.TiltFac.x;
				}
			}
			Delta.z = D1.Atoms[i].AtomCoord.z - D1.Atoms[j].AtomCoord.z;
			Delta.z = PBCDistance(Delta.z,Lz);
			Distance = DistanceSquare(Delta.x,Delta.y,Delta.z);

			if (sqrt(Distance) < L_half)
				{
					BinIndex = int(sqrt(Distance)/BinSize);
					rdf[BinIndex] = rdf[BinIndex] + 2.0;

				}		

		}
	}
	//Normalizing RDF
	double N_Ideal,Shell_Volume,r;
	double V1,V2,delV;
	double Volume = Lx * Ly*Lz; // 2d VOLUME
	Rho = (D1.NumberofAtoms)/Volume;

	for (int BinIndex = 0; BinIndex< NumberofBins;BinIndex++)
	{
		double r1 = BinIndex*BinSize;
		double r2 = r1 + BinSize;
		delV = 3.14159*((r1+r2)*(r1+r2))*BinSize;
		rdf[BinIndex] /= (Rho*D1.NumberofAtoms*delV);

	}

}