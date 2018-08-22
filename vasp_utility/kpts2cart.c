#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[]) {

	FILE *kpoints, *outcar;
	FILE *kpts2cart;

	int i, j;
	char buffer[512];

	int nbands, nkpts;
	double **kpts, b[3][3];
	char str[14] = "energy-cutoff";


        if(argc != 2) {
                printf("Usage: kpts2cart <# of k-points>\n");
                return 0;
        }
        else {
		if((kpoints = fopen("KPOINTS", "r")) == NULL) {
			printf("KPOINTS is not found.\n");
			return 0;
		}
		if((outcar = fopen("OUTCAR", "r")) == NULL) {
			printf("OUTCAR is not found.\n");
			return 0;
		}
		nkpts = atoi(argv[1]);
	}

	fscanf(outcar, "%s", buffer);
	while(memcmp(buffer, str, strlen(str)) != 0) fscanf(outcar, "%s", buffer);
	for(i=0; i<3; i++) fgets(buffer, 512, outcar);
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) fscanf(outcar, "%s", buffer);
		for(j=0; j<3; j++) fscanf(outcar, "%lf", &b[i][j]);
	}

	fclose(outcar);


	kpts2cart = fopen("KPOINTS_CART", "w");

	for(i=0; i<2; i++) {
		fgets(buffer, 512, kpoints);
		fputs(buffer, kpts2cart);
	}
	fgets(buffer, 512, kpoints);
	fprintf(kpts2cart, "Cartesian\n");

	kpts = (double**)malloc(sizeof(double*)*nkpts);
	for(i=0; i<nkpts; i++) kpts[i] = (double*)malloc(sizeof(double)*4);

	for(i=0; i<nkpts; i++) {
		for(j=0; j<4; j++) fscanf(kpoints, "%lf", &kpts[i][j]);
	}
	for(i=0; i<nkpts; i++) {
		fprintf(kpts2cart, "%14.10f %14.10f %14.10f %14.10f\n", kpts[i][0]*b[0][0]+kpts[i][1]*b[1][0]+kpts[i][2]*b[2][0], kpts[i][0]*b[0][1]+kpts[i][1]*b[1][1]+kpts[i][2]*b[2][1], kpts[i][0]*b[0][2]+kpts[i][1]*b[1][2]+kpts[i][2]*b[2][2], kpts[i][3]);
        }

	free(kpts);
	fclose(kpoints);
	fclose(kpts2cart);

	return 0;

}
