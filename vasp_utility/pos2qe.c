#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {

	FILE *poscar, *qe;

	int i, j;
	char buffer[512];

	int pc, nti, tni, *ni;
	char **na;
	double sf, a[3][3], ***coord;


	if(argc != 2) {
		printf("Usage: pos2qe <# of types of ions>\n");
		return 0;
	}
	else {
		if((poscar = fopen("POSCAR", "r")) == NULL) {
			printf("POSCAR is not found.\n");
			return 0;
		}
		nti = atoi(argv[1]);
	}

	ni = (int*)malloc(sizeof(int)*nti);

	na = (char**)malloc(sizeof(char*)*nti);
	for(i=0; i<nti; i++) na[i] = (char*)malloc(sizeof(char)*3);

	fgets(buffer, 512, poscar);
	fscanf(poscar, "%lf", &sf);
	for(i=0; i<3; i++) for(j=0; j<3; j++) fscanf(poscar, "%lf", &a[i][j]);
	for(i=0; i<nti; i++) fscanf(poscar, "%s", na[i]);
	for(i=0; i<nti; i++) fscanf(poscar, "%d", &ni[i]);
	tni = 0;
	for(i=0; i<nti; i++) tni += ni[i];
	fscanf(poscar, "%s", buffer);
	if(buffer[0] == 'C') pc = 0;
	else if(buffer[0] == 'D') pc = 1;

	coord = (double***)malloc(sizeof(double**)*nti);
	for(i=0; i<nti; i++) coord[i] = (double**)malloc(sizeof(double*)*tni);
	for(i=0; i<nti; i++) for(j=0; j<tni; j++) coord[i][j] = (double*)malloc(sizeof(double)*3);

	for(i=0; i<nti; i++) {
		for(j=0; j<ni[i]; j++) fscanf(poscar, "%lf %lf %lf", &coord[i][j][0], &coord[i][j][1], &coord[i][j][2]);
	}

	fclose(poscar);


	qe = fopen("geometry.dat", "w");

	fprintf(qe, " CELL_PARAMETERS (angstrom)\n");
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) fprintf(qe, "  %13.9f", sf*a[i][j]);
		fprintf(qe, "\n");
	}
	if(pc == 0) fprintf(qe, " ATOMIC_POSITIONS (angstrom)\n");
	if(pc == 1) fprintf(qe, " ATOMIC_POSITIONS (crystal)\n");
	for(i=0; i<nti; i++) {
                for(j=0; j<ni[i]; j++) {
			fprintf(qe, "   %s  %13.9f  %13.9f  %13.9f\n", na[i], coord[i][j][0], coord[i][j][1], coord[i][j][2]);
		}
	}

	free(ni);
	free(na);
	free(coord);
	fclose(qe);
	
	return 0;

}
