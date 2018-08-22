/* 
  For Brillouin zone sampling (n x n square shape)
  co[0] Gamma point
  co[1] kx point
  co[2] ky point
*/

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[]) {

	FILE *kpts_out;

	int i, j, k, div, qf;
	double co[3][2]={0}, **co1, ***co2;

	if (argc != 5) {
		printf("Usage: kpts_gen_2d <kx> <ky> <# of division> <quarter(0)/full(1)>\n");
		return 1;
	}
	else {
		kpts_out = fopen("kpts.out", "w");
		co[1][0] = atof(argv[1]);
		co[2][1] = atof(argv[2]);
		div = atoi(argv[3]);
		qf = atoi(argv[4]);
	}
	co1 = (double**)malloc(sizeof(double*)*div);
	for(i=0; i<div; i++) co1[i] = (double*)malloc(sizeof(double)*2);
	co2 = (double***)malloc(sizeof(double**)*div);
	for(i=0; i<div; i++) co2[i] = (double**)malloc(sizeof(double*)*div);
	for(i=0; i<div; i++) for(j=0; j<div; j++) co2[i][j] = (double*)malloc(sizeof(double)*2);

	for(i=0; i<div; i++) {
		for(j=0; j<2; j++) {
			if(i == 0) co1[i][j] = co[i][j];
			else co1[i][j] = i*((co[1][j]-co[0][j])/(div-1)) + co[0][j];
		}
	}

	for(i=0; i<div; i++) {
		for(j=0; j<div; j++) {
			for(k=0; k<2; k++) {
				if(j == 0) co2[i][j][k] = co1[i][k];
				else co2[i][j][k] = j*((co[2][k]-co[0][k])/(div-1)) + co1[i][k];
			}
		}
	}

        for(i=0; i<div; i++) {
                for(j=0; j<div; j++) {
			fprintf(kpts_out, "  %14.10f  %14.10f  %14.10f  %14.10f\n", co2[i][j][0], co2[i][j][1], 0.0, 0.0);
		}
	}
	if(qf == 1) {
		for(i=0; i<div; i++) {
			for(j=0; j<div; j++) {
				fprintf(kpts_out,"  %14.10f  %14.10f  %14.10f  %14.10f\n", -co2[i][j][0], co2[i][j][1], 0.0, 0.0);
			}
		}
		for(i=0; i<div; i++) {
			for(j=0; j<div; j++) {
				fprintf(kpts_out,"  %14.10f  %14.10f  %14.10f  %14.10f\n", co2[i][j][0], -co2[i][j][1], 0.0, 0.0);
			}
		}
		for(i=0; i<div; i++) {
			for(j=0; j<div; j++) {
				fprintf(kpts_out,"  %14.10f  %14.10f  %14.10f  %14.10f\n", -co2[i][j][0], -co2[i][j][1], 0.0, 0.0);
			}
		}
	}
	free(co1);
	free(co2);
	fclose(kpts_out);
	
	return 0;
}
