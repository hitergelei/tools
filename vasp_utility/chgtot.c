#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

	FILE *parchg, *chgtot;

	int  i, j, k, count;
	char buffer[512];

	int ntyps, nions, ngrid[3];
	double co, chg, *sum;


	if((parchg = fopen("PARCHG", "r")) == NULL) {
		printf("PARCHG is not found.\n");
		return 0;
	}

	for(i=0; i<2; i++) fgets(buffer, 512, parchg);
	for(i=0; i<1; i++) fscanf(parchg, "%lf", &co);
	for(i=0; i<4; i++) fgets(buffer, 512, parchg);
	count = 0;
	while(1) {
		fgets(buffer, 512, parchg);
		if(buffer[0]=='D' || buffer[0]=='C') break;
		count++;
	}
	if(count == 0) count = 7;
	else if(count == 1) count = 8;

	rewind(parchg);

	ntyps = 0;
	for(i=0; i<count-2; i++) fgets(buffer, 512, parchg);
	while(1) {
		fscanf(parchg, "%s", buffer);
		if(buffer[0]=='C' || buffer[0]=='D') break;
		ntyps++;
	}

	rewind(parchg);

	nions = 0;
	for(i=0; i<count-2; i++) fgets(buffer, 512, parchg);
	for(i=0; i<ntyps; i++) {
		fscanf(parchg, "%d", &j);
		nions += j;
	}
	for(i=0; i<nions+3; i++) fgets(buffer, 512, parchg);
	for(i=0; i<3; i++) fscanf(parchg, "%d", &ngrid[i]);


	chgtot = fopen("pchg.dat", "w");

	sum = (double*)malloc(sizeof(double)*ngrid[0]);

	sum[0] = 0.0;
	for(i=0; i<ngrid[1]; i++) {
		for(j=0; j<ngrid[0]; j++) {
			fscanf(parchg, "%lf", &chg);
			sum[j] += chg;
		}
	}
	for(i=0; i<ngrid[2]-1; i++) {
		for(j=0; j<ngrid[1]; j++) {
			for(k=0; k<ngrid[0]; k++) {
				fscanf(parchg, "%lf", &chg);
				sum[k] += chg;
			}
		}
	}
	for(i=0; i<ngrid[0]; i++) fprintf(chgtot, "  %14.10f  %14.10f\n", co/ngrid[0]*(i+1), sum[i]/(ngrid[1]*ngrid[2]));

	free(sum);
        fclose(parchg);
        fclose(chgtot);

}
