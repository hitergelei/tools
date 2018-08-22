#include <stdio.h>
#include <stdlib.h>

#define VACCUM 15

int main(int argc, char* argv[]) {

	FILE *poscar, *slab;

	int i, j, k;
	char buffer[512];

	int nti, tni, axis, ntimes, *ni;
	char **na;
	double sf, a[3][3];
	double ***coord;


	if (argc != 4) {
		printf("Usage: pos2slab <# of types of ions> <axis(a,b,c = 1,2,3)> <n times>\n");
		return 0;
	}
	else {
		if((poscar = fopen("POSCAR", "r")) == NULL) {
			printf("POSCAR in not found.\n");
			return 0;
		}
		nti = atoi(argv[1]);
		axis = atoi(argv[2]);
		ntimes = atoi(argv[3]);
        }

	slab = fopen("POSCAR_SLAB", "w");

	fgets(buffer, 512, poscar);
	fputs(buffer, slab);
	fscanf(poscar, "%lf", &sf);
	fprintf(slab, "  %.10f\n", sf);
	for(i=0; i<3; i++) for(j=0; j<3; j++) fscanf(poscar, "%lf", &a[i][j]);
	if(axis == 1) {
		fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", ntimes*a[0][0]+VACCUM, a[0][1], a[0][2]);
		fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", a[1][0], a[1][1], a[1][2]);
		fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", a[2][0], a[2][1], a[2][2]);
	}
	else if (axis == 2) {
		fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", a[0][0], a[0][1], a[0][2]);
		fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", a[1][0], ntimes*a[1][1]+VACCUM, a[1][2]);
		fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", a[2][0], a[2][1], a[2][2]);
	}
	else if (axis == 3) {
		fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", a[0][0], a[0][1], a[0][2]);
		fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", a[1][0], a[1][1], a[1][2]);
		fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", a[2][0], a[2][1], ntimes*a[2][2]+VACCUM);
	}

	na = (char**)malloc(sizeof(char*)*nti);
	for(i=0; i<nti; i++) na[i] = (char*)malloc(sizeof(char)*3);

	ni = (int*)malloc(sizeof(int)*nti);

	for(i=0; i<nti; i++) fscanf(poscar, "%s", na[i]);
	for(i=0; i<nti; i++) fprintf(slab, "  %s", na[i]);
	fprintf(slab, "\n");
	for(i=0; i<nti; i++) fscanf(poscar, "%d", &ni[i]);
	for(i=0; i<nti; i++) fprintf(slab, "  %d", ntimes*ni[i]);
	tni = 0;
	for(i=0; i<nti; i++) tni += ni[i];
	fprintf(slab, "\n");
	fscanf(poscar, "%s", buffer);
	fprintf(slab, "%s\n", buffer);

	coord = (double***)malloc(sizeof(double**)*nti);
	for(i=0; i<nti; i++) coord[i] = (double**)malloc(sizeof(double*)*tni);
	for(i=0; i<nti; i++) for(j=0; j<tni; j++) coord[i][j] = (double*)malloc(sizeof(double)*3);

	for(i=0; i<nti; i++) 
		for(j=0; j<ni[i]; j++) fscanf(poscar, "%lf %lf %lf", &coord[i][j][0], &coord[i][j][1], &coord[i][j][2]);
	if(axis == 1) {
		for(i=0; i<nti; i++) {
			for(k=0; k<ntimes; k++) {
		                for(j=0; j<ni[i]; j++) {
	        	                fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", k*a[0][0]+coord[i][j][0], coord[i][j][1], coord[i][j][2]);
				}
			}
		}
	}
	if(axis == 2) {
		for(i=0; i<nti; i++) {
			for(k=0; k<ntimes; k++) {
				for(j=0; j<ni[i]; j++) {
					fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", coord[i][j][0], k*a[1][1]+coord[i][j][1], coord[i][j][2]);
				}
			}
		}
	}
	if(axis == 3) {
		for(i=0; i<nti; i++) {
			for(k=0; k<ntimes; k++) {
				for(j=0; j<ni[i]; j++) {
					fprintf(slab, "  %14.10f  %14.10f  %14.10f\n", coord[i][j][0], coord[i][j][1], k*a[2][2]+coord[i][j][2]);
				}
			}
		}
	}

	free(na);
	free(ni);
	free(coord);
	fclose(slab);
	fclose(poscar);
	
	return 0;

}
