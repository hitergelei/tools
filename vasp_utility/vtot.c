#include <stdio.h>

int main(int argc, char *argv[]) {

	FILE *locpot, *vtot;

	int  i, j, count;
	char buffer[512];

	int nti, ni, ngrid[3];
	double z_coord, potential, sum;


	if((locpot = fopen("LOCPOT", "r")) == NULL) {
		printf("LOCPOT is not found.\n");
		return 0;
	}

	for(i=0; i<4; i++) fgets(buffer, 512, locpot);  
	for(i=0; i<3; i++) fscanf(locpot, "%lf", &z_coord);
	for(i=0; i<2; i++) fgets(buffer, 512, locpot);
	count = 0;
	while(1) {
		fgets(buffer, 512, locpot);
		if(buffer[0]=='D' || buffer[0]=='C') break;
		count++;
	}
	if(count == 0) count = 7;
	else if(count == 1) count = 8;

	rewind(locpot);

	for(i=0; i<count-2; i++) fgets(buffer, 512, locpot);
	nti = 0;
	while(1) {
		fscanf(locpot, "%s", buffer);
		if(buffer[0]=='C' || buffer[0]=='D') break;
		nti++;
	}

	rewind(locpot);

	ni = 0;
	for(i=0; i<count-2; i++) fgets(buffer, 512, locpot);
	for(i=0; i<nti; i++) {
		fscanf(locpot, "%d", &j);
		ni += j;
	}
	for(i=0; i<ni+3; i++) fgets(buffer, 512, locpot);
	for(i=0; i<3; i++) fscanf(locpot, "%d", &ngrid[i]);


	vtot = fopen("vtot.dat", "w");

	for(i=0; i<ngrid[2]; i++) {
		sum = 0;
		for(j=0; j<ngrid[0]*ngrid[1]; j++) {
			fscanf(locpot, "%lf", &potential);
			sum += potential;
		}
		fprintf(vtot, "  %14.10f  %14.10f\n", z_coord/ngrid[2]*(i+1), sum/(ngrid[0]*ngrid[1]));
	}

	fclose(locpot);
	fclose(vtot);

	return 0;

}
