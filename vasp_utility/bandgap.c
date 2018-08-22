#include <stdio.h>
#include <string.h>

typedef struct {
	int nkpts;
	int nbands;
	int ispin;
} PARA;

typedef struct {
	double kpts[2][3];
	double bands[2];
} KE;

int ReadOUTCAR (char *str);
void KEp (FILE *outcar, int ispin);
KE KEdata (FILE *outcar, PARA p);
void Result (KE dat);

int main(int argc, char* argv[]) {

	FILE *outcar;

	int i;
	char buffer[512];

	PARA p;
	KE dat;


	p.nkpts = ReadOUTCAR("NKPTS");
	p.nbands = ReadOUTCAR("NBANDS=");
	p.ispin = ReadOUTCAR("ISPIN");

	outcar = fopen("OUTCAR", "r");

	KEp(outcar, p.ispin);

	for(i=0; i<3; i++) fscanf(outcar, "%s", buffer);
	dat = KEdata(outcar, p);
	if(p.ispin == 2) printf(" ------  Majority spin bands  ------\n");
	Result(dat);
	if(p.ispin == 2) {
		for(i=0; i<6; i++) fscanf(outcar, "%s", buffer);
		dat = KEdata(outcar, p);
		printf(" ------  Minority spin bands  ------\n");
		Result(dat);
	}

	fclose(outcar);

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
	else fscanf(fp, "%s %d", buffer, &out);

	fclose(fp);

	return out;

}

void KEp (FILE *outcar, int ispin) {

	int i;
	char buffer[512];

	char str[10]="E-fermi";


	while(memcmp(buffer, str, strlen(str)) != 0) fscanf(outcar, "%s", buffer);
	if(ispin == 1) {
		for(i=0; i<3; i++) fgets(buffer, 512, outcar);
	}
	else if(ispin == 2) {
		for(i=0; i<5; i++) fgets(buffer, 512, outcar);
	}

}

KE KEdata (FILE *outcar, PARA p) {

	int i, j, k, count;
	char buffer[512];

	KE out;
	double ktmp[3], btmp[3];


	for(i=0; i<3; i++) fscanf(outcar, "%lf", &ktmp[i]);
	for(i=0; i<2; i++) for(j=0; j<3; j++) out.kpts[i][j] = ktmp[j];
	for(i=0; i<2; i++) fgets(buffer, 512, outcar);
	count = 0;
	for(i=0; i<p.nbands; i++) {
		for(j=0; j<3; j++) fscanf(outcar, "%lf", &btmp[j]);
		if(btmp[2] >= 0.5) {
			out.bands[0] = btmp[1];
		} // VBM of first k-point
		else if (count == 0) {
			out.bands[1] = btmp[1];
			++count;
		} // CBM of first k-point
	}
	for(i=0; i<p.nkpts-1; i++) {
		for(j=0; j<3; j++) fscanf(outcar, "%s", buffer);
		for(j=0; j<3; j++) fscanf(outcar, "%lf", &ktmp[j]);
		for(j=0; j<2; j++) fgets(buffer, 512, outcar);
		for(j=0; j<p.nbands; j++) {
			for(k=0; k<3; k++) fscanf(outcar, "%lf", &btmp[k]);
			if((btmp[2] >= 0.5) && (btmp[1] > out.bands[0])) {
				out.bands[0] = btmp[1];
				for(k=0;k<3;k++) out.kpts[0][k] = ktmp[k];
			}
			if((btmp[2] < 0.5) && (btmp[1] < out.bands[1])) {
				out.bands[1] = btmp[1];
				for(k=0; k<3; k++) out.kpts[1][k] = ktmp[k];
			}
		}
	}

	return out;

}

void Result (KE dat) {

	char buffer[10];


	printf(" CBM is at %8.6f %8.6f %8.6f\n", dat.kpts[1][0], dat.kpts[1][1], dat.kpts[1][2]);
	printf(" VBM is at %8.6f %8.6f %8.6f\n", dat.kpts[0][0], dat.kpts[0][1], dat.kpts[0][2]);
	if(dat.kpts[0][0] == dat.kpts[1][0] && dat.kpts[0][1] == dat.kpts[1][1] && dat.kpts[0][2] == dat.kpts[1][2]) {
		 strcpy(buffer, "Direct");
	}
	else strcpy(buffer, "Indirect");
	printf(" The center of the bandgap : %8.6f\n", (dat.bands[0]+dat.bands[1])/2);
	printf(" BANDGAP(%s) : %8.6f\n", buffer, dat.bands[1]-dat.bands[0]);

}
