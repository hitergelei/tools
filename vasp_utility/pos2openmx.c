#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {

	FILE *poscar, *openmx;

	int i, j, k;
	char buffer[512];

	int nti, tni, *ni;
	char **na;
	double sf, a[3][3];
	double **ne, ***coord;


	if (argc != 2) {
                printf("Usage: pos2openmx <# of types of ions>\n");
                return 0;
        }
        else {
		if((poscar = fopen("POSCAR", "r")) == NULL) {
			printf("POSCAR is not found.\n");
			return 0;
		}
		nti = atoi(argv[1]);
        }

	ne = (double**)malloc(sizeof(double*)*nti);
	for(i=0; i<nti; i++) ne[i] = (double*)malloc(sizeof(double)*2);

	na = (char**)malloc(sizeof(char*)*nti);
	for(i=0; i<nti; i++) na[i] = (char*)malloc(sizeof(char)*3);

	ni = (int*)malloc(sizeof(int)*nti);

	printf("Enter the number of electrons of each atoms for spin up & down\n");
	printf("(e.g.) 1.0 0.0 for H atom\n");
	for(i=0; i<nti; i++) scanf("%lf %lf", &ne[i][0], &ne[i][1]);

	fgets(buffer, 512, poscar);
	fscanf(poscar, "%lf", &sf);
	for(i=0; i<3; i++) for(j=0; j<3; j++) fscanf(poscar, "%lf", &a[i][j]);
	for(i=0; i<nti; i++) fscanf(poscar, "%s", na[i]);
	for(i=0; i<nti; i++) fscanf(poscar, "%d", &ni[i]);
	tni = 0;
	for(i=0; i<nti; i++) tni += ni[i];

	coord = (double***)malloc(sizeof(double**)*nti);
	for(i=0; i<nti; i++) coord[i] = (double**)malloc(sizeof(double*)*tni);
	for(i=0; i<nti; i++) for(j=0; j<tni; j++) coord[i][j] = (double*)malloc(sizeof(double)*3);

	fscanf(poscar, "%s", buffer);
	for(i=0; i<nti; i++) {
		for(j=0; j<ni[i]; j++) fscanf(poscar, "%lf %lf %lf", &coord[i][j][0], &coord[i][j][1], &coord[i][j][2]);
	}

	fclose(poscar);


	openmx = fopen("geometry.dat", "w");

	fprintf(openmx, "Atoms.Number  %d\n", tni);
	fprintf(openmx, "Atoms.SpeciesAndCoordinates.Unit  Ang  # Ang|AU\n");
	fprintf(openmx, "<Atoms.SpeciesAndCoordinates\n");
	k = 1;
	for(i=0; i<nti; i++) {
                for(j=0; j<ni[i]; j++) {
			fprintf(openmx,"  %3.0d  %s  %14.10f  %14.10f  %14.10f  %1.1f  %1.1f\n", k, na[i], coord[i][j][0], coord[i][j][1], coord[i][j][2], ne[i][0], ne[i][1]);
			k++;
		}
	}
	fprintf(openmx, "Atoms.SpeciesAndCoordinates>\n");
	fprintf(openmx, "Atoms.UnitVectors.Unit            Ang  # Ang|AU\n");
	fprintf(openmx, "<Atoms.UnitVectors\n");
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) fprintf(openmx, "  %14.10f", sf*a[i][j]);
		fprintf(openmx, "\n");
	}
	fprintf(openmx, "Atoms.UnitVectors>\n");

	free(ne);
	free(na);
	free(ni);
	free(coord);
	fclose(openmx);
	
	return 0;

}
