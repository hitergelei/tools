#include <stdio.h>
#include <stdlib.h>

int main(void) {

	FILE *contcar, *poscar;

	int i, j;
	char buffer[512];

	int ena;
	char cod;

	int nti, tni, *ni;
	char **na;
	double sf, a[3][3], ***coord;


	if((contcar = fopen("CONTCAR", "r")) == NULL) {
		printf("CONTCAR is not found.\n");
		return 0;
	}
	for(i=0; i<8; i++) {
		fgets(buffer, 512, contcar);
		if(i>=6 && (buffer[0]=='C' || buffer[0]=='D')) {
			ena = i + 1;
			cod = buffer[0];
			break;
		}
	}

	rewind(contcar);

	if(ena == 7) {
		for(i=0; i<5; i++) fgets(buffer, 512, contcar);
	}
	else if(ena == 8) {
		for(i=0; i<6; i++) fgets(buffer, 512, contcar);
	}
	for(i=0;  ; i++) {
		fscanf(contcar, "%s", buffer);
		if(buffer[0]=='C' || buffer[0]=='D') break;
	}
	nti = i;

	rewind(contcar);

	poscar = fopen("POSCAR", "w");

	fgets(buffer, 512, contcar);
	fputs(buffer, poscar); // NAME OF THE SYSTEM
	fscanf(contcar, "%lf", &sf);
	fprintf(poscar, "  %.10f\n", sf); // SCALING FACTOR
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) fscanf(contcar, "%lf", &a[i][j]);
	}
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) fprintf(poscar, "  %20.16f", a[i][j]);
		fprintf(poscar, "\n");
	} // LATTICE VECTORS

	ni = (int*)malloc(sizeof(int)*nti);

	na = (char**)malloc(sizeof(char*)*nti);
	for(i=0; i<nti; i++) na[i] = (char*)malloc(sizeof(char)*3);

	if(ena == 7) {
		for(i=0; i<nti; i++) fscanf(contcar, "%d", &ni[i]);
		for(i=0; i<nti; i++) fprintf(poscar, "  %d", ni[i]);
		fprintf(poscar, "\n");
	} // NUMBER OF IONS (ABSENCE of NAME OF ATOMS)
	else if(ena == 8) {
		for(i=0; i<nti; i++) fscanf(contcar, "%s", na[i]);
		for(i=0; i<nti; i++) fscanf(contcar, "%d", &ni[i]);
		for(i=0; i<nti; i++) fprintf(poscar, "  %s", na[i]);
		fprintf(poscar, "\n");
		for(i=0; i<nti; i++) fprintf(poscar, "  %d", ni[i]);
		fprintf(poscar, "\n");
		
	} // NAME OF ATOMS & NUMBER OF IONS
	tni = 0;
	for(i=0; i<nti; i++) tni += ni[i]; // TOTAL NUMBER OF IONS

	coord = (double***)malloc(sizeof(double**)*nti);
	for(i=0; i<nti; i++) coord[i] = (double**)malloc(sizeof(double*)*tni);
	for(i=0; i<nti; i++) for(j=0; j<tni; j++) coord[i][j] = (double*)malloc(sizeof(double)*3);

	fscanf(contcar, "%s", buffer);
	for(i=0; i<nti; i++)
		for(j=0; j<ni[i]; j++) fscanf(contcar, "%lf %lf %lf", &coord[i][j][0], &coord[i][j][1], &coord[i][j][2]);
	if(cod == 'D') {
		fprintf(poscar, "Cartesian\n");
		for(i=0; i<nti; i++)
			for(j=0; j<ni[i]; j++)
				fprintf(poscar, "  %20.16f  %20.16f  %20.16f\n", sf*(coord[i][j][0]*a[0][0]+coord[i][j][1]*a[1][0]+coord[i][j][2]*a[2][0]), sf*(coord[i][j][0]*a[0][1]+coord[i][j][1]*a[1][1]+coord[i][j][2]*a[2][1]), sf*(coord[i][j][0]*a[0][2]+coord[i][j][1]*a[1][2]+coord[i][j][2]*a[2][2]));
	}
	if(cod == 'C') {
		fprintf(poscar, "Cartesian\n");
		for(i=0; i<nti; i++)
			for(j=0; j<ni[i]; j++)
				fprintf(poscar, "  %20.16f  %20.16f  %20.16f\n", coord[i][j][0], coord[i][j][1], coord[i][j][2]);
	}

	free(ni);
	free(na);
	free(coord);
	fclose(poscar);
	fclose(contcar);
	
	return 0;

}
