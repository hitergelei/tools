#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[]) {

	FILE *pos, *xsf;
	
	int i=0, j=0, k=0;
	char buffer[512];

	int nti, tni, ni[10];
	char **na, **na_tmp;
	double a[3][3], **coord;

	
	if (argc != 2) {
		printf("Usage: xsf2pos <the name of the system>\n");
		return 0;
	}
        else {
		if((xsf = fopen("POSCAR.xsf", "r")) != NULL) {}
		else if((xsf = fopen("CONTCAR.xsf", "r")) != NULL) {}
		else {
			printf(".xsf file is not found.\n");
			return 0;
		}
	}

	for(i=0; i<2; i++) fgets(buffer, 512, xsf);
	for(i=0; i<3; i++) for(j=0; j<3; j++) fscanf(xsf, "%lf", &a[i][j]); // READ LATTICE CONSTANTS
	for(i=0; i<2; i++) fgets(buffer, 512, xsf);
	fscanf(xsf, "%d %s", &tni, buffer); // READ TOTAL NUMBER OF IONS

	na = (char**)malloc(sizeof(char*)*tni);
	for(i=0; i<tni; i++) na[i] = (char*)malloc(sizeof(char)*3);

	na_tmp = (char**)malloc(sizeof(char*)*tni);
	for(i=0; i<tni; i++) na_tmp[i] = (char*)malloc(sizeof(char)*3);

	coord = (double**)malloc(sizeof(double*)*tni);
	for(i=0; i<tni; i++) coord[i] = (double*)malloc(sizeof(double)*3);

	for(i=0; i<tni; i++) {
		fscanf(xsf, "%s %lf %lf %lf", na_tmp[i], &coord[i][0], &coord[i][1], &coord[i][2]);
		k++;
		if(i!=0 && (strcmp(na_tmp[i-1], na_tmp[i])!=0)) {
			strcpy(na[j], na_tmp[i-1]);
			ni[j] = k - 1;
			k = 1;
			j++;
		}
		if(i == tni-1) {
			ni[j] = k;
			strcpy(na[j], na_tmp[i]);
		}
	} // READ NAME & COORDINATES OF EACH ION
	nti = j + 1; // NUMBER OF TYPES OF IONS

	free(na_tmp);
	fclose(xsf);


	pos = fopen("POSCAR", "w");

	fprintf(pos, "%s\n", argv[1]);
	fprintf(pos, "  1.\n");
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) fprintf(pos, "  %14.10f", a[i][j]);
		fprintf(pos, "\n");
	}
	for(i=0; i<nti; i++) fprintf(pos, "  %s", na[i]);
	fprintf(pos, "\n");
	for(i=0; i<nti; i++) fprintf(pos, "  %d", ni[i]);
	fprintf(pos, "\n");
	fprintf(pos, "Cartesian\n");
	for(i=0; i<tni; i++) fprintf(pos, "  %14.10f  %14.10f  %14.10f\n", coord[i][0], coord[i][1], coord[i][2]);

	free(na);
	free(coord);
	fclose(pos);

	return 0;

}
