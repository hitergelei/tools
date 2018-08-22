/*
 This program makes 2D banddata for Mathematica.
*/
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {

	FILE *eigenval, *outcar;
	FILE *banddat;

	int i, j;
	char buffer[512];

	int nbands, nkpts, ispin;
	double Ef, Emin, Emax;
	double b[3][3], **kpts, ***Eband;


	if(argc != 4) {
		printf("Usage: 2dbd <E_min> <E_max> <Fermi energy>\n");
		return 0;
	}
	else {
		if((eigenval = fopen("EIGENVAL", "r")) == NULL) {
			printf("EIGENVAL is not found.\n");
			return 0;
		}
		if((outcar = fopen("OUTCAR", "r")) == NULL) {
			printf("OUTCAR is not found.\n");
			return 0;
		}
		Emin = atof(argv[1]);
		Emax = atof(argv[2]);
		Ef = atof(argv[3]);
	}

	for(i=0; i<4; i++) fscanf(eigenval, "%d", &ispin);
	for(i=0; i<5; i++) fgets(buffer, 512, eigenval);
	fscanf(eigenval, "%s %d %d", buffer, &nkpts, &nbands);
	for(i=0; i<2; i++) fgets(buffer, 512, eigenval);

	kpts = (double**)malloc(sizeof(double*)*nkpts);
	for(i=0; i<nkpts; i++) kpts[i] = (double*)malloc(sizeof(double)*3);

	Eband = (double***)malloc(sizeof(double**)*nkpts);
	for(i=0; i<nkpts; i++) Eband[i] = (double**)malloc(sizeof(double*)*nbands);
	for(i=0; i<nkpts; i++) for(j=0; j<nbands; j++) Eband[i][j] = (double*)malloc(sizeof(double)*2);

	for(i=0; i<nkpts; i++) {
		fscanf(eigenval, "%lf %lf %lf %s", &kpts[i][0], &kpts[i][1], &kpts[i][2], buffer);
		if(ispin == 2)
			for(j=0; j<nbands; j++) fscanf(eigenval, "%s %lf %lf", buffer, &Eband[i][j][0], &Eband[i][j][1]);
		else if(ispin == 1)
			for(j=0; j<nbands; j++) fscanf(eigenval, "%s %lf", buffer, &Eband[i][j][0]);
		for(j=0; j<2; j++) fgets(buffer, 512, eigenval);
	}

	fclose(eigenval);

	while(1) {
		fscanf(outcar, "%s", buffer);
		if(buffer[5]=='y' && buffer[6]=='-' && buffer[7]=='c') {
			for(i=0; i<3; i++) fgets(buffer, 512, outcar);
			break;
		}
	}
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) fscanf(outcar, "%s", buffer);
		for(j=0; j<3; j++) fscanf(outcar, "%lf", &b[i][j]);
	}

	fclose(outcar);

	if(ispin == 1) {
		banddat = fopen("2d.bd", "w");

		for(i=0; i<nkpts; i++) {
			for(j=0; j<nbands; j++) {
				if(Eband[i][j][0]-Ef>Emin && Eband[i][j][0]-Ef<Emax)
					fprintf(banddat, "%f %f      %f\n", kpts[i][0]*b[0][0]+kpts[i][1]*b[1][0]+kpts[i][2]*b[2][0], kpts[i][0]*b[0][1]+kpts[i][1]*b[1][1]+kpts[i][2]*b[2][1], Eband[i][j][0]-Ef);
			}
		}

		fclose(banddat);
	}

	if(ispin == 2) {
		banddat = fopen("2dup.bd", "w");

		for(i=0; i<nkpts; i++) {
			for(j=0; j<nbands; j++) {
				if(Eband[i][j][0]-Ef>Emin && Eband[i][j][0]-Ef<Emax)
					fprintf(banddat, "%f	%f	%f\n", kpts[i][0]*b[0][0]+kpts[i][1]*b[1][0]+kpts[i][2]*b[2][0], kpts[i][0]*b[0][1]+kpts[i][1]*b[1][1]+kpts[i][2]*b[2][1], Eband[i][j][0]-Ef);
			}
		}

		fclose(banddat);

		banddat = fopen("2ddw.bd", "w");

		for(i=0; i<nkpts; i++) {
			for(j=0; j<nbands; j++) {
				if(Eband[i][j][0]-Ef>Emin && Eband[i][j][0]-Ef<Emax)
					fprintf(banddat, "%f	%f	%f\n", kpts[i][0]*b[0][0]+kpts[i][1]*b[1][0]+kpts[i][2]*b[2][0], kpts[i][0]*b[0][1]+kpts[i][1]*b[1][1]+kpts[i][2]*b[2][1], Eband[i][j][1]-Ef);
			}
		}

		fclose(banddat);
	}

	free(kpts);
	free(Eband);

	return 0;

}
