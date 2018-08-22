#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
	int mode;
	int *range;
	int nedos;
	int ispin;
	int soc;
	int lorbit;
	int nit;
} PARA;

int ReadOUTCAR (char *str);
void WriteSumDOS (int argc, PARA p);

int main(int argc, char* argv[]) {

	int i;
	char buffer[512];

	PARA p;
	double **pdos, **sum;


	if(argc < 4) {
		printf("Usage: sumdos [0|1] [0:n1 n2]|[1:n1 n2 n3 ...]\n");
		printf("(e.g.) sumdos 0 1 10 (from dos1 to dos10) || sumdos 1 1 2 3 4 (dos1, dos2, dos3, and dos4)\n");
		return 0;
	}
	
	p.mode = atoi(argv[1]);
	if(p.mode == 0) {
		p.range = (int*)malloc(sizeof(int)*2);
		p.range[0] = atoi(argv[2]);
		p.range[1] = atoi(argv[3]);
	} // MODE 0
	else if(p.mode == 1) {
		p.range = (int*)malloc(sizeof(int)*(argc-2));
		for(i=0; i<argc-2; i++) p.range[i] = atoi(argv[i+2]);
	} // MODE 1
	else {
		printf("Choose mode 0 or 1.\n");
		return 0;
	}

	p.nedos = ReadOUTCAR("NEDOS");
	p.ispin = ReadOUTCAR("ISPIN");
	p.soc = ReadOUTCAR("LNONCOLLINEAR");
	p.lorbit = ReadOUTCAR("LORBIT");

	if(p.lorbit==10 && p.ispin==1 && p.soc==0) p.nit = 4;
	else if(p.lorbit==10 && p.ispin==2 && p.soc==0) p.nit = 7;
	else if(p.lorbit==10 && p.soc==1) p.nit = 13;
	else if(p.lorbit==11 && p.ispin==1 && p.soc==0) p.nit = 10;
	else if(p.lorbit==11 && p.ispin==2 && p.soc==0) p.nit = 19;
	else if(p.lorbit==11 && p.soc==1) p.nit = 37;

	WriteSumDOS(argc, p);

	free(p.range);

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

void WriteSumDOS (int argc, PARA p) {

	FILE *dos, *sumdos;

	int i, j, k;
	char buffer[512];

	int ini, fin;
	double **pdos, **sum;


	pdos = (double**)malloc(sizeof(double*)*p.nedos);
	for(i=0; i<p.nedos; i++) pdos[i] = (double*)malloc(sizeof(double)*p.nit);

	sum = (double**)malloc(sizeof(double*)*p.nedos);
	for(i=0; i<p.nedos; i++) sum[i] = (double*)malloc(sizeof(double)*p.nit);

	if(p.mode == 0) {
		ini = p.range[0];
		fin = p.range[1] + 1;
	}
	else if(p.mode == 1) {
		ini = 0;
		fin = argc - 2;
	}
	for(i=ini; i<fin; i++) {
		if(p.mode == 0) sprintf(buffer, "dos%d", i);
		else if(p.mode == 1) sprintf(buffer, "dos%d", p.range[i]);

		if((dos = fopen(buffer, "r")) == NULL) {
			if(p.mode == 0) printf("dos%d is not found.\n", i);
			else if(p.mode == 1) printf("dos%d is not found.\n", p.range[i]);
		}

		for(j=0; j<p.nedos; j++) {
			for(k=0; k<p.nit; k++) fscanf(dos, "%lf", &pdos[j][k]);
			for(k=1; k<p.nit; k++) sum[j][k] += pdos[j][k];
		}

       		fclose(dos);
	}

	sumdos = fopen("sumdos", "w");

	for(i=0; i<p.nedos; i++) {
		fprintf(sumdos, "  %14.10f", pdos[i][0]);
		for(j=1; j<p.nit; j++) fprintf(sumdos, "  %14.10f", sum[i][j]);
		fprintf(sumdos, "\n");
	}

	free(pdos);
	free(sum);
	fclose(sumdos);

}
