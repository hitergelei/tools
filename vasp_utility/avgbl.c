#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double BondLength (double **co, int i, int j);

int main(int argc, char* argv[]) {

	FILE *poscar;

	int i, j, k;
	char buffer[512];

	int nt, *na;
	double cutrad, **co, *leng;
	double sum;

	
	if(argc < 3) {
		printf("Usage: avgbl <# of types of atoms> <cutoff radius>\n");
		return 0;
	}
	else {
		nt = atoi(argv[1]);
		cutrad = atof(argv[2]);
	}

	if((poscar = fopen("POSCAR", "r")) == NULL) {
		printf("POSCAR is not found.\n");
		return 0;
	}

	for(i=0; i<6; i++) fgets(buffer, 512, poscar);

	na = (int*)malloc(sizeof(int)*nt);

	for(i=0; i<nt; i++) fscanf(poscar, "%d", &na[i]);
	for(i=1; i<nt; i++) na[0] += na[i];
	for(i=0; i<2; i++) fgets(buffer, 512, poscar);

	co = (double**)malloc(sizeof(double*)*na[0]);
	for(i=0; i<na[0]; i++) co[i] = (double*)malloc(sizeof(double)*3);

	for(i=0; i<na[0]; i++) {
		for(j=0; j<3; j++) fscanf(poscar, "%lf", &co[i][j]);
	}

	fclose(poscar);

	leng = (double*)malloc(sizeof(double)*(na[0]*na[0]));

	k=0;
	for(i=0; i<na[0]-1; i++) {
		for(j=0; j<na[0]; j++) {
			if(j > i) {
				leng[k] = BondLength(co, i, j);
				if(leng[k] > cutrad) --k;
				++k;
			}
		}
	}

	sum=0;
	for(i=0; i<k; i++) sum += leng[i];
	printf(" Average bond length: %f A\n", sum/(double)k);

	free(na);
	free(co);
	free(leng);

	return 0;

}

double BondLength (double **co, int i, int j) {
	
	double bl;


	bl = sqrt(pow((co[i][0]-co[j][0]), 2.0) + pow((co[i][1]-co[j][1]), 2.0) + pow((co[i][2]-co[j][2]), 2.0));
	
	return bl;

}
