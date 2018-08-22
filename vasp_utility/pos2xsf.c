#include <stdio.h>
#include <stdlib.h>

int main(void) {

	FILE *poscar, *xsf;

	int i, j;
	char buffer[512];

	int ena;
	char cod;

	int nti, tni, *ni;
	char **na; 
	double sf, a[3][3], ***coord;


	if((poscar = fopen("POSCAR", "r")) == NULL) {
		printf("POSCAR is not found.\n");
		return 0;
	}
	for(i=0; i<8; i++) {
		fgets(buffer, 512, poscar);
		if(i>=6 && (buffer[0]=='C' || buffer[0]=='D')) {
			ena = i + 1;
			cod = buffer[0];
			break;
		}
	}

	rewind(poscar);

	if(ena == 7) for(i=0; i<5; i++) fgets(buffer, 512, poscar);
	else if(ena == 8) for(i=0; i<6; i++) fgets(buffer, 512, poscar);
	for(i=0;  ; i++) {
		fscanf(poscar, "%s", buffer);
		if(buffer[0]=='C' || buffer[0]=='D') break;
	}
	nti = i;

	rewind(poscar);

	fgets(buffer, 512, poscar);
	fscanf(poscar, "%lf", &sf);
	for(i=0; i<3; i++) for(j=0; j<3; j++) fscanf(poscar, "%lf", &a[i][j]);

	ni = (int*)malloc(sizeof(int)*nti);

	na = (char**)malloc(sizeof(char*)*nti);
	for(i=0; i<nti; i++) na[i] = (char*)malloc(sizeof(char)*3);

	if(ena == 7) for(i=0; i<nti; i++) fscanf(poscar, "%d", &ni[i]);
	else if(ena == 8) {
		for(i=0; i<nti; i++) fscanf(poscar, "%s", na[i]);
		for(i=0; i<nti; i++) fscanf(poscar, "%d", &ni[i]);
	}

	coord = (double***)malloc(sizeof(double**)*nti);
	for(i=0; i<nti; i++) coord[i] = (double**)malloc(sizeof(double*)*ni[i]);
	for(i=0; i<nti; i++) for(j=0; j<ni[i]; j++) coord[i][j] = (double*)malloc(sizeof(double)*3);

	for(i=0; i<2; i++) fgets(buffer, 512, poscar);
	for(i=0; i<nti; i++) {
		for(j=0; j<ni[i]; j++) fscanf(poscar,"%lf %lf %lf", &coord[i][j][0], &coord[i][j][1], &coord[i][j][2]);
	}

	fclose(poscar);


	xsf = fopen("POSCAR.xsf", "w");

	fprintf(xsf, "CRYSTAL\nPRIMVEC\n");
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) fprintf(xsf, "  %14.10f", sf*a[i][j]);
		fprintf(xsf, "\n");
	}
	fprintf(xsf, "PRIMCOORD\n");
	tni = 0;
	for(i=0; i<nti; i++) tni += ni[i];
	fprintf(xsf, "   %d  1\n", tni);
	for(i=0; i<nti; i++) {
		for(j=0; j<ni[i]; j++) {
			fprintf(xsf,"  %3s  %14.10f  %14.10f  %14.10f\n", na[i], coord[i][j][0], coord[i][j][1], coord[i][j][2]);
		}
	}

	free(ni);
	free(na);
	free(coord);
	fclose(xsf);
	
	return 0;

}
