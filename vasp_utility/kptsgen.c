/* 
  For 0-weight(fake) SC procedure (band structure with hybrid fuctional)
  After running this program, type "cat kpts.out > KPOINTS"
*/
#include <stdio.h>
#include <stdlib.h>

int main(void) {

	FILE *kpts_in, *kpts_out;

	int i, j, k;
	int nlines, ndiv;

	double **coord1, ***coord2;


	if((kpts_in = fopen("kpts.in", "r")) == NULL) {
		printf("kpts.in is not found. The contents of kpts.in are as follows.\n");
		printf("<# of lines> <# of division of line>\n");
		printf("The list of k-points\n");
		return 0;
	}

	fscanf(kpts_in, "%d %d", &nlines, &ndiv);

	coord1 = (double**)malloc(sizeof(double*)*(nlines+1));
	for(i=0; i<nlines+1; i++) coord1[i] = (double*)malloc(sizeof(double)*3);

	for(i=0; i<nlines+1; i++) fscanf(kpts_in, "%lf %lf %lf", &coord1[i][0], &coord1[i][1], &coord1[i][2]);

	fclose(kpts_in);

	coord2 = (double***)malloc(sizeof(double**)*nlines);
	for(i=0; i<nlines; i++) coord2[i] = (double**)malloc(sizeof(double*)*(nlines*ndiv));
	for(i=0; i<nlines; i++) for(j=0; j<nlines*ndiv; j++) coord2[i][j] = (double*)malloc(sizeof(double)*3);

	for(i=0; i<3; i++) coord2[0][0][i] = coord1[0][i];
	for(i=0; i<nlines; i++) {
		for(j=0; j<ndiv; j++) {
			for(k=0; k<3; k++) {
				if(i==0 && j==0) coord2[i][j][k] += (coord1[i+1][k]-coord1[i][k])/(double)ndiv;
				else {
					if(j == 0) coord2[i][j][k] = coord2[i-1][ndiv-1][k] + (coord1[i+1][k]-coord1[i][k])/(double)ndiv;
					else {
						coord2[i][j][k] = coord2[i][j-1][k] + (coord1[i+1][k]-coord1[i][k])/(double)ndiv;
						if(coord2[i][j][k]<1E-12 && coord2[i][j][k]>-1E-12) coord2[i][j][k] = 0.0;
					}
				}
			}
		}
	}


	kpts_out = fopen("kpts.out", "w");

	fprintf(kpts_out,"  %14.10f  %14.10f  %14.10f  %14.10f\n", coord1[0][0], coord1[0][1], coord1[0][2], 0.0);
	for(i=0; i<nlines; i++) { 
		for(j=0; j<ndiv; j++) {
			fprintf(kpts_out,"  %14.10f  %14.10f  %14.10f  %14.10f\n", coord2[i][j][0], coord2[i][j][1], coord2[i][j][2], 0.0);
		}
	}

	fclose(kpts_out);
	free(coord1);
	free(coord2);

	return 0;

}
