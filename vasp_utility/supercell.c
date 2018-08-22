#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {

	FILE *poscar, *sposcar;

	int i, j, k, l, m;
	char buffer[512];

	int ntimes;
	char sob;

	int ena;
	char cod;

	int nti, *ni;
	char **na, ns[256];
	double sf, a[3][3], ***coord;


	if (argc != 3) {
		printf("Usage: supercell <n times> <slab(y) or bulk(n)>\n");
		return 0;
	}
	else {
		if((poscar = fopen("POSCAR", "r")) == NULL) {
			printf("POSCAR is not found.\n");
			return 0;
		}
		ntimes = atoi(argv[1]);
		sob = argv[2][0];
	}

	for(i=0; i<8; i++) {
		fgets(buffer, 512, poscar);
		if(i>=6 && (buffer[0]=='C' || buffer[0]=='D')) {
			ena = i + 1;
			cod = buffer[0];
			break;
		}
	}
	if(cod == 'D') {
		printf("Change diret coordinates to cartesian coordinates!\n");
		return 0;
	}

	rewind(poscar);

	if(ena == 7) for(i=0; i<5; i++) fgets(buffer, 512, poscar);
	else if(ena == 8) for(i=0; i<6; i++) fgets(buffer, 512, poscar);
	for(i=0;  ; i++) {
		fscanf(poscar, "%s", buffer);
		if(buffer[0] == 'C') break;
	}
	nti = i;

	rewind(poscar);

	ni = (int*)malloc(sizeof(int)*nti);

	na = (char**)malloc(sizeof(char*)*nti);
	for(i=0; i<nti; i++) na[i] = (char*)malloc(sizeof(char)*3);

	fgets(ns, 256, poscar);
	fscanf(poscar, "%lf", &sf);
	for(i=0; i<3; i++) for(j=0; j<3; j++) fscanf(poscar, "%lf", &a[i][j]);
	if(ena == 8) for(i=0; i<nti; i++) fscanf(poscar, "%s", na[i]);
	for(i=0; i<nti; i++) fscanf(poscar, "%d", &ni[i]);
//	tni = 0;
//	for(i=0; i<nti; i++) tni += ni[i];

	coord = (double***)malloc(sizeof(double**)*nti);
	for(i=0; i<nti; i++) coord[i] = (double**)malloc(sizeof(double*)*ni[i]);
	for(i=0; i<nti; i++) for(j=0; j<ni[i]; j++) coord[i][j] = (double*)malloc(sizeof(double)*3);

	fscanf(poscar, "%s", buffer);
	for(i=0; i<nti; i++) {
		for(j=0; j<ni[i]; j++) fscanf(poscar, "%lf %lf %lf", &coord[i][j][0], &coord[i][j][1], &coord[i][j][2]);
	}

	fclose(poscar);


	sposcar = fopen("SPOSCAR", "w");

	/********************slab = no********************/
	if(sob == 'n') {
		fprintf(sposcar, "%s", ns);
		fprintf(sposcar, "  %.10f\n", sf);
		fprintf(sposcar, "  %14.10f  %14.10f  %14.10f\n", ntimes*a[0][0], ntimes*a[0][1], ntimes*a[0][2]);
     		fprintf(sposcar, "  %14.10f  %14.10f  %14.10f\n", ntimes*a[1][0], ntimes*a[1][1], ntimes*a[1][2]);
		fprintf(sposcar, "  %14.10f  %14.10f  %14.10f\n", ntimes*a[2][0], ntimes*a[2][1], ntimes*a[2][2]);
		if(ena == 8) {
			for(i=0; i<nti; i++) fprintf(sposcar, "  %s", na[i]);
			fprintf(sposcar, "\n");
		}
		for(i=0; i<nti; i++) fprintf(sposcar, "  %d", ntimes*ntimes*ntimes*ni[i]);
		fprintf(sposcar, "\n");
		fprintf(sposcar, "Cartesian\n");
		for(i=0; i<nti; i++) 
			for(j=0; j<ni[i]; j++)
				for(k=0; k<ntimes; k++)
					for(l=0; l<ntimes; l++)
						for(m=0; m<ntimes; m++) fprintf(sposcar, "  %14.10f  %14.10f  %14.10f\n", k*a[0][0]+l*a[1][0]+m*a[2][0]+coord[i][j][0], k*a[0][1]+l*a[1][1]+m*a[2][1]+coord[i][j][1], k*a[0][2]+l*a[1][2]+m*a[2][2]+coord[i][j][2]);
	}

	/********************slab = yes********************/
	else {
		fprintf(sposcar, "%s", ns);
		fprintf(sposcar, "  %.10f\n", sf);
		fprintf(sposcar, "  %14.10f  %14.10f  %14.10f\n", ntimes*a[0][0], ntimes*a[0][1], ntimes*a[0][2]);
		fprintf(sposcar, "  %14.10f  %14.10f  %14.10f\n", ntimes*a[1][0], ntimes*a[1][1], ntimes*a[1][2]);
		fprintf(sposcar, "  %14.10f  %14.10f  %14.10f\n", a[2][0], a[2][1], a[2][2]);
		if(ena == 8) {
			for(i=0; i<nti; i++) fprintf(sposcar, "  %s", na[i]);
			fprintf(sposcar, "\n");
		}
		for(i=0; i<nti; i++) fprintf(sposcar, "  %d", ntimes*ntimes*ni[i]);
		fprintf(sposcar, "\n");
		fprintf(sposcar, "Cartesian\n");
		for(i=0; i<nti; i++) 
			for(j=0; j<ni[i]; j++)
				for(k=0; k<ntimes; k++)
					for(l=0; l<ntimes; l++) fprintf(sposcar, "  %14.10f  %14.10f  %14.10f\n", k*a[0][0]+l*a[1][0]+coord[i][j][0], k*a[0][1]+l*a[1][1]+coord[i][j][1], k*a[0][2]+l*a[1][2]+coord[i][j][2]);
	}

	free(ni);
	free(na);
	free(coord);
	fclose(sposcar);
	
	return 0;

}
