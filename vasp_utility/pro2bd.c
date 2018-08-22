#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
        double Ef;
        int nkpts;
        int nbands;
        int nlines;
	int nions;
        int soc;
	int lorbit;
	int nsi;
	int *lsi;
} PARA;

typedef struct {
        int hf;
        int nnzw;
} NZW;

typedef struct {
	double Eband;
	double kptsco;
} BS;

int ReadOUTCAR (char *str);
double*** ReadKPTS (int nlines, int nkpts, NZW n);
double* CalcDKPTS (int nlines, double ***kpts);
double*** MakeArry (int s1, int s2, int s3);
void ReadPROCAR(FILE *procar, FILE *pro2bd, PARA p, BS b);
NZW NumZeroWeight (int nkpts);

int main(int argc, char* argv[]) {

        FILE *procar;
	FILE *pro2bd;

	int i, j, k, m;
	char buffer[512];

	PARA p;
	NZW n;
	BS b;
	double ***kpts, *dkpts;


	if (argc != 5) {
		printf("Usage: pro2bd <selected ions> <Fermi energy> <# of lines>\n");
		return 0;
	}
	else {
		if((procar = fopen("PROCAR", "r")) == NULL) {
			printf("PROCAR is not found.\n");
			return 0;
		}
		if(*argv[1] == *argv[2]) p.nsi = 1;
		else p.nsi = atoi(argv[2]) - atoi(argv[1]) + 1;

		p.lsi = (int*)malloc(sizeof(int)*p.nsi);

		p.lsi[0] = atoi(argv[1]);
		for(i=0; i<p.nsi-1; i++) p.lsi[i+1] = p.lsi[i] + 1;
		p.Ef = atof(argv[3]);
		p.nlines = atoi(argv[4]);
	}

/** READ OUTCAR **/
	p.nkpts = ReadOUTCAR("NKPTS");
	p.nbands = ReadOUTCAR("NBANDS=");
	p.nions = ReadOUTCAR("NIONS");
	p.soc = ReadOUTCAR("LSORBIT");
	p.lorbit = ReadOUTCAR("LORBIT");
	n = NumZeroWeight(p.nkpts);
	if(n.hf == 1) p.nkpts -= n.nnzw;
	if(n.hf == 0) kpts = MakeArry(p.nlines, (int)p.nkpts/p.nlines, 4);
	if(n.hf == 1) kpts = MakeArry(p.nlines, (int)(p.nkpts-1)/p.nlines+1, 4);
	kpts = ReadKPTS(p.nlines, p.nkpts, n);
/** CLOSE OUTCAR **/

	dkpts = (double*)malloc(sizeof(double)*(p.nlines+1));

	dkpts = CalcDKPTS(p.nlines, kpts);

	free(kpts);

/** READ PROCAR & WRITE proj.bd **/
	m = 1;
	b.kptsco = 0.0;

	pro2bd = fopen("proj.bd", "w");

	for(i=0; i<3; i++) fgets(buffer, 512, procar);
	if(n.hf == 1) {
		for(i=0; i<n.nnzw+1; i++) {
			fscanf(procar, "%s", buffer);
			while(memcmp(buffer, "k-point", strlen("k-point")) != 0) fscanf(procar, "%s", buffer);
		}
	}
	for(i=0; i<p.nkpts; i++) {
		fgets(buffer, 512, procar);
		for(j=0; j<p.nbands; j++) {
			for(k=0; k<4; k++) fscanf(procar, "%s", buffer);
			fscanf(procar, "%lf", &b.Eband);
			fgets(buffer, 512, procar);
			if(p.lorbit == 10) for(k=0; k<5; k++) fscanf(procar, "%s", buffer);
			else if(p.lorbit == 11) for(k=0; k<11; k++) fscanf(procar, "%s", buffer);
			ReadPROCAR(procar, pro2bd, p, b);
			if(p.soc == 0) fgets(buffer, 512, procar);
			else if(p.soc == 1) for(k=0; k<3*(p.nions+1)+1; k++) fgets(buffer, 512, procar);
		}
		b.kptsco += dkpts[0];
		if(n.hf == 0) {
			if(i == m*(p.nkpts/p.nlines)-1) {
				b.kptsco -= dkpts[0];
				for(j=0; j<p.nlines-1; j++) dkpts[j] = dkpts[j+1];
				m++;
			}
		}
		if(n.hf == 1) {
			if(i == m*((p.nkpts-1)/p.nlines)-1) {
				for(j=0; j<p.nlines-1; j++) dkpts[j] = dkpts[j+1];
				m++;
			}
		}
		fprintf(pro2bd, "\n");
		for(j=0; j<2; j++) fgets(buffer, 512, procar);
	}

	free(p.lsi);
	free(dkpts);
	fclose(pro2bd);
	fclose(procar);
/** CLOSE PROCAR & proj.bd **/

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
	if(strcmp(str, "NBANDS=") == 0) fscanf(fp, "%d", &out);
	else if(strcmp(str, "LSORBIT") == 0) {
		for(i=0; i<2; i++) fscanf(fp, "%s", buffer);
		if(strcmp(buffer, "T") == 0) out = 1;
		else out = 0;
	}
	else fscanf(fp, "%s %d", buffer, &out);

	fclose(fp);

	return out;

}

double*** ReadKPTS (int nlines, int nkpts, NZW n) {

	FILE *fp;

	int i, j, k;
	char buffer[512];

	double ***out;


	if((fp = fopen("OUTCAR", "r")) == NULL) {
		printf("OUTCAR is not found.\n");
		return 0;
	}
	fscanf(fp, "%s", buffer);
	while(memcmp(buffer, "2pi/SCALE", strlen("2pi/SCALE")) != 0) fscanf(fp, "%s", buffer);
	if(n.hf == 0) {
		fgets(buffer, 512, fp);
		out = MakeArry(nlines, (int)nkpts/nlines, 4);
		for(i=0; i<nlines; i++) {
			for(j=0; j<nkpts/nlines; j++) {
				for(k=0; k<4; k++) fscanf(fp, "%lf", &out[i][j][k]);
			}
		}
	}
	else if(n.hf == 1) {
		for(i=0; i<n.nnzw+1; i++) fgets(buffer, 512, fp);
		out = MakeArry(nlines, (int)(nkpts-1/nlines)+1, 4);
		for(i=0; i<nlines; i++) {
			for(j=0; j<(nkpts-1)/nlines; j++) {
				for(k=0; k<4; k++) fscanf(fp, "%lf", &out[i][j][k]);
				if((i == nlines-1) && (j == (nkpts-1)/nlines-1)) {
					for(k=0; k<4; k++) fscanf(fp, "%lf", &out[i][j+1][k]);
				}
			}
		}
	}

	fclose(fp);

	return out;

}

double* CalcDKPTS (int nlines, double ***kpts) {

	int i;
	double *out;


	out = (double*)malloc(sizeof(double)*(nlines+1));

	for(i=0; i<nlines; i++) out[i] = sqrt((pow((kpts[i][0][0]-kpts[i][1][0]), 2.0) + pow((kpts[i][0][1]-kpts[i][1][1]), 2.0) + pow((kpts[i][0][2]-kpts[i][1][2]), 2.0)));

	return out;

}

double*** MakeArry (int s1, int s2, int s3) {

	int i, j;
	double ***out;


	out = (double***)malloc(sizeof(double**)*s1);
	for(i=0; i<s1; i++) out[i] = (double**)malloc(sizeof(double*)*s2);
	for(i=0; i<s1; i++) for(j=0; j<s2; j++) out[i][j] = (double*)malloc(sizeof(double)*s3);

	return out;

}

void ReadPROCAR (FILE *procar, FILE *pro2bd, PARA p, BS b) {

	int i, j, k;
	char buffer[512];

	int num;
	double *sum, **tmp;


	if(p.lorbit == 10) num = 4;
	else if(p.lorbit == 11) num = 10;

	sum = (double*)malloc(sizeof(double)*(num+1));

	tmp = (double**)malloc(sizeof(double*)*(p.nions+2));
	for(i=0; i<p.nions+1; i++) tmp[i] = (double*)malloc(sizeof(double)*(num+1));

	for(i=0; i<p.nions+1; i++) {
		fscanf(procar, "%s", buffer);
		for(j=0; j<num; j++) fscanf(procar, "%lf", &tmp[i][j]);
	}
	for(i=0; i<num; i++) sum[i] = 0;
	for(i=0; i<p.nions; i++) {
		for(j=0; j<p.nsi; j++) {
			if(i == p.lsi[j]-1) {
				for(k=0; k<num; k++) sum[k] += tmp[i][k];
			}
		}
		if(i == p.nions-1) {
			fprintf(pro2bd, "%12.8f  %12.8f", b.kptsco, b.Eband-p.Ef);
			for(j=0; j<num; j++) fprintf(pro2bd, "  %12.8f", sum[j]);
			fprintf(pro2bd, "\n");
		}
	}

	free(sum);
	free(tmp);

}

NZW NumZeroWeight (int nkpts) {

	FILE *fp;

	int i, j;
	char buffer[512];
	double tmp[4];

	NZW out;


	if((fp = fopen("OUTCAR", "r")) == NULL) {
		printf("OUTCAR is not found.\n");
	}

	fscanf(fp, "%s", buffer);
	while(memcmp(buffer, "2pi/SCALE", strlen("2pi/SCALE")) != 0) fscanf(fp, "%s", buffer);
	fgets(buffer, 512, fp);
	out.nnzw = 0;
	for(i=0; i<nkpts; i++) {
		for(j=0; j<4; j++) fscanf(fp, "%lf", &tmp[j]);
		if(tmp[3] == 0) {
			out.hf = 1;
			break;
		}
		else {
			out.nnzw++;
			out.hf = 0;
		}
	}

	fclose(fp);

	return out;

}
