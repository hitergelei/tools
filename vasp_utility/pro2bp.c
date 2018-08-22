#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>

#define NUM_FILES 10

typedef struct {
	int nkpts;
	int nlines;
	double Emin;
	double Emax;
} PARA;

typedef struct {
	int hf;
	int nnzw;
} NZW;

int ReadOUTCAR (char *str);
double*** ReadKPTS (int nlines, int nkpts, NZW n);
NZW NumZeroWeight (int nkpts);
double* CalcDKPTS (int nlines, double ***kpts);
double*** MakeArry (int s1, int s2, int s3);

int main(int argc, char* argv[]) {

	DIR *wd;
	struct dirent *entry;

	FILE *pro2bp, *tic;

	int i, j, k;

	PARA p;
	NZW n;
	int fleng, nfiles;
	double ***kpts, *dkpts, *hsptsco;
	char hspts[20], fnames[NUM_FILES][128]={0};


	if (argc != 4) {
		printf("Usage: pro2bp <E_min> <E_max> <high-symmetry points>\n");
		return 0;
	}
	else {
		p.Emin = atof(argv[1]); 
		p.Emax = atof(argv[2]);
		strcpy(hspts, argv[3]);
		p.nlines = strlen(hspts) - 1;
        }

	p.nkpts = ReadOUTCAR("NKPTS");
	n = NumZeroWeight(p.nkpts);
	if(n.hf == 0) kpts = MakeArry(p.nlines, (int)p.nkpts/p.nlines, 4);
	if(n.hf == 1) {
		p.nkpts -= n.nnzw;
		kpts = MakeArry(p.nlines, (int)(p.nkpts-1)/p.nlines+1, 4);
	}
	kpts = ReadKPTS(p.nlines, p.nkpts, n);

	dkpts = (double*)malloc(sizeof(double)*(p.nlines+1));

	dkpts = CalcDKPTS(p.nlines, kpts);

	free(kpts);

	hsptsco = (double*)malloc(sizeof(double)*(p.nlines+1));

	hsptsco[0] = 0.0;
	if(n.hf == 0) for(i=0; i<p.nlines; i++) hsptsco[i+1] = hsptsco[i] + dkpts[i]*(p.nkpts/p.nlines-1);
	if(n.hf == 1) for(i=0; i<p.nlines; i++) hsptsco[i+1] = hsptsco[i] + dkpts[i]*((p.nkpts-1)/p.nlines);

/** WRITE Projbp.in & tic.dat **/
	if(p.nlines != 1) tic = fopen("tic.dat", "w");
	pro2bp = fopen("Projbp.in", "w");

	fprintf(pro2bp, "set size 0.5,0.5\n");
	fprintf(pro2bp, "set terminal postscript enhanced 12\n");
	fprintf(pro2bp, "set output 'proj_band.eps'\n");
	fprintf(pro2bp, "set nokey\n");
	fprintf(pro2bp, "set ylabel 'Energy (eV)'\n");
	fprintf(pro2bp, "set xrange [0:%f]\n", hsptsco[p.nlines]);
	fprintf(pro2bp, "set yrange [%f:%f]\n", p.Emin, p.Emax);
	fprintf(pro2bp, "set ytics %f,%f,%f\n", p.Emin, (p.Emax-p.Emin)/2, p.Emax);
	fprintf(pro2bp, "set pointsize 0.2\n");
	fprintf(pro2bp, "set bmargin 3\n");
	fprintf(pro2bp, "set xzeroaxis\n");
	for(i=0; i<p.nlines+1; i++) {
		if(hspts[i]=='G' && i==0)
			fprintf(pro2bp, "set xtics (""\"{/Symbol %c}\""" %f,", hspts[i], hsptsco[i]);
		if(hspts[i]!='G' && i==0)
			fprintf(pro2bp, "set xtics (""\"%c\""" %f,", hspts[i], hsptsco[i]);
		if(hspts[i]=='G' && i!=0 && i!=p.nlines)
			fprintf(pro2bp, " ""\"{/Symbol %c}\""" %f,", hspts[i], hsptsco[i]);
		if(hspts[i]!='G' && i!=0 && i!=p.nlines)
			fprintf(pro2bp, " ""\"%c\""" %f,", hspts[i], hsptsco[i]);
		if(hspts[i]=='G' && i==p.nlines)
			fprintf(pro2bp, " ""\"{/Symbol %c}\""" %f)\n", hspts[i], hsptsco[i]);
		if(hspts[i]!='G' && i==p.nlines)
			fprintf(pro2bp, " ""\"%c\""" %f)\n", hspts[i], hsptsco[i]);
	}

	wd = opendir("./");

	entry = readdir(wd);
	for(i=0;  ;  ) {
		while(1) {
			j = 0;
			fleng = 0;
			while(entry->d_name[j] != '\0') {
				fleng++;
				j++;
			}
			if(fleng >= 4) {
				if(entry->d_name[fleng-1]=='d' && entry->d_name[fleng-2]=='b' && entry->d_name[fleng-3]=='.') {
	                                for(k=0;k<fleng;k++) fnames[i][k] = entry->d_name[k];
					i++;
		                }
			}
			entry = readdir(wd);
			break;
		}
		if(entry == NULL) break;
	}
	nfiles = i;

	closedir(wd);


	if(nfiles != 0) fprintf(pro2bp, "plot	\'%s\' using 1:2:($3*0.01) with circles lc rgb ""\"black\""" ,\\\n", fnames[0]);
	for(i=0; i<nfiles-1; i++) {
		fprintf(pro2bp, "	\'%s\' using 1:2:($3*0.01) with circles lc rgb ""\"black\""" ,\\\n", fnames[i+1]);
	}
	if(p.nlines != 1) {
		fprintf(pro2bp, "	'tic.dat' using 1:2 with lines\n");
		for(i=0; i<p.nlines-1; i++) for(j=-50; j<51; j++) fprintf(tic, "%f  %d\n", hsptsco[i+1], j);
	}

	free(dkpts);
	free(hsptsco);
        fclose(pro2bp);
	fclose(tic);

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

NZW NumZeroWeight (int nkpts) {

	FILE *fp;

	int i, j;
	char buffer[512];

	NZW out;
	double tmp[4];


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
