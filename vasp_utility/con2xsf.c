#include <stdio.h>
#include <stdlib.h>

int main(void) {

	FILE *contcar, *xsf;

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
			ena = i+1;
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

	fgets(buffer, 512, contcar);
	fscanf(contcar, "%lf", &sf);
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) fscanf(contcar, "%lf", &a[i][j]);
	}

	ni = (int*)malloc(sizeof(int)*nti);

	na = (char**)malloc(sizeof(char*)*nti);
	for(i=0; i<nti; i++) na[i] = (char*)malloc(sizeof(char)*3);

	if(ena == 7) {
		for(i=0; i<nti; i++) fscanf(contcar, "%d", &ni[i]);
	}
        else if(ena == 8) {
                for(i=0; i<nti; i++) fscanf(contcar, "%s", na[i]);
                for(i=0; i<nti; i++) fscanf(contcar, "%d", &ni[i]);
        }
	tni = 0;
	for(i=0; i<nti; i++) tni += ni[i];
	for(i=0; i<2; i++) fgets(buffer, 512, contcar);

	coord = (double***)malloc(sizeof(double**)*nti);
	for(i=0; i<nti; i++) coord[i] = (double**)malloc(sizeof(double*)*tni);
	for(i=0; i<nti; i++) for(j=0; j<tni; j++) coord[i][j] = (double*)malloc(sizeof(double)*3);

	for(i=0; i<nti; i++) { 
		for(j=0; j<ni[i]; j++) fscanf(contcar, "%lf %lf %lf", &coord[i][j][0], &coord[i][j][1], &coord[i][j][2]);
	}

	fclose(contcar);


	xsf = fopen("CONTCAR.xsf", "w");

	fprintf(xsf, "CRYSTAL\nPRIMVEC\n");
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) fprintf(xsf, "  %14.10f", sf*a[i][j]);
		fprintf(xsf, "\n");
	}
	fprintf(xsf, "PRIMCOORD\n");
	fprintf(xsf, "   %d  1\n", tni);
	for(i=0; i<nti; i++) {
		for(j=0; j<ni[i]; j++)
			fprintf(xsf, "  %3s  %14.10f  %14.10f  %14.10f\n", na[i], sf*(coord[i][j][0]*a[0][0]+coord[i][j][1]*a[1][0]+coord[i][j][2]*a[2][0]), sf*(coord[i][j][0]*a[0][1]+coord[i][j][1]*a[1][1]+coord[i][j][2]*a[2][1]), sf*(coord[i][j][0]*a[0][2]+coord[i][j][1]*a[1][2]+coord[i][j][2]*a[2][2]));
	}

	free(ni);
	free(na);
	free(coord);
	fclose(xsf);
	
	return 0;

}
