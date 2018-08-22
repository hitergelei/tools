#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int ReadOUTCAR (char *str);

int main(int argc, char* argv[]) {

	FILE *doscar, *dos;

	int i, j, k;
	char buffer[512];

	int nedos, ispin, soc, lorbit, ni, nit;
	double Ef, *tdos, *pdos;


	if((doscar = fopen("DOSCAR", "r")) == NULL) {
		printf("DOSCAR is not found.\n");
		return 0;
	}

	ispin = ReadOUTCAR("ISPIN");
	soc = ReadOUTCAR("LNONCOLLINEAR");
	lorbit = ReadOUTCAR("LORBIT");

	fscanf(doscar, "%d", &ni);
	for(i=0; i<5; i++) fgets(buffer, 512, doscar);
	fscanf(doscar, "%s %s %d %lf %s", buffer, buffer, &nedos, &Ef, buffer);


	dos = fopen("dos0", "w");

	tdos = (double*)malloc(sizeof(double)*(2*ispin+1));

	for(i=0; i<nedos; i++) {
		for(j=0; j<(2*ispin+1); j++) fscanf(doscar, "%lf", &tdos[j]);
		fprintf(dos, " %14.10f", tdos[0]-Ef);
		for(j=1; j<(2*ispin+1); j++) fprintf(dos, "  %14.10f", tdos[j]);
		fprintf(dos, "\n");
	}

	free(tdos);
	fclose(dos);

	if(lorbit==10 && ispin==1 && soc==0) nit = 4;
	else if(lorbit==10 && ispin==2 && soc==0) nit = 7;
	else if(lorbit==10 && soc==1) nit = 13;
	else if(lorbit==11 && ispin==1 && soc==0) nit = 10;
	else if(lorbit==11 && ispin==2 && soc==0) nit = 19;
	else if(lorbit==11 && soc==1) nit = 37;

	pdos = (double*)malloc(sizeof(double)*nit);

	for(i=0; i<ni; i++) {
		sprintf(buffer, "dos%d", i+1);
		dos = fopen(buffer, "w");

		for(j=0; j<2; j++) fgets(buffer, 512, doscar);
		for(j=0; j<nedos; j++) {
			for(k=0; k<nit; k++) fscanf(doscar, "%lf", &pdos[k]);
			fprintf(dos, "  %14.10f", pdos[0]-Ef);
			for(k=1; k<nit; k++) fprintf(dos, "  %14.10f", pdos[k]);
			fprintf(dos, "\n");
		}

		fclose(dos);
	}

	free(pdos);
	fclose(doscar);

	return 0;

}

int ReadOUTCAR (char *str) {

	FILE *fp;

	int i;
	char buffer[512];

	int out;


	if((fp = fopen("OUTCAR", "r")) == NULL) {
		printf("OUTCAR is not found.\n");
		return 0;
	}

	fscanf(fp, "%s", buffer);
	while(memcmp(buffer, str, strlen(str)) != 0) fscanf(fp, "%s", buffer);
	if(strcmp(str, "LNONCOLLINEAR") == 0) {
		for(i=0; i<2; i++) fscanf(fp, "%s", buffer);
		if(strcmp(buffer, "T") == 0) out = 1;
		else out = 0;
	}
	else fscanf(fp, "%s %d", buffer, &out);

	fclose(fp);

	return out;

}
