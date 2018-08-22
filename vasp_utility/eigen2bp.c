#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
	double Emin;
	double Emax;
	double Ef;
	int nkpts;
	int nbands;
	int nlines;
	int ispin;
	int soc;
	char hspts[20];
} PARA;

typedef struct {
	int hf;
	int nnzw;
} NZW;

int ReadOUTCAR (char *str);
double*** ReadKPTS (int nlines, int nkpts, NZW n);
double*** ReadEIGENVAL (PARA p, NZW n);
double*** MakeArry (int s1, int s2, int s3);
double* CalcDKPTS (int nlines, double ***kpts);
void WriteBanddat (PARA p, double *kptsco, double ***Eband);
void WriteGNUscript (PARA p, double *hsptsco);
NZW NumZeroWeight (int nkpts);

int main(int argc, char* argv[]) {

	int i, j, m;

	PARA p;
	NZW n;
	double *dkpts, *hsptsco, *kptsco;
	double ***kpts, ***Eband;


        if (argc != 5) {
		printf("Usage: eigen2bp <E_min> <E_max> <high-symmetry points> <Fermi energy>\n");
		printf("(e.g.) eigen2bp -1 1 MKGM 1.00\n");
		return 0;
        }
        else {
                p.Emin = atof(argv[1]);
                p.Emax = atof(argv[2]);
		strcpy(p.hspts, argv[3]);
                p.nlines = strlen(p.hspts)-1;
	 	p.Ef = atof(argv[4]);
        }
	
/** READ OUTCAR **/
	p.nkpts = ReadOUTCAR("NKPTS");
	p.nbands = ReadOUTCAR("NBANDS=");
	p.ispin = ReadOUTCAR("ISPIN");
	p.soc = ReadOUTCAR("LSORBIT");
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

	kptsco = (double*)malloc(sizeof(double)*(p.nkpts+1));

	kptsco[0] = 0.0;
	m = 1;
        for(i=1; i<p.nkpts; i++) {
                kptsco[i] = kptsco[i-1] + dkpts[0];
		if(n.hf == 0) {
	                if(i == m*(p.nkpts/p.nlines)) {
				kptsco[i] -= dkpts[0];
	                        for(j=0; j<p.nlines-1; j++) dkpts[j] = dkpts[j+1];
	                        m++;
			}
		}
		else if(n.hf == 1) {
			if(i == m*((p.nkpts-1)/p.nlines)) {
				for(j=0; j<p.nlines-1; j++) dkpts[j] = dkpts[j+1];
				m++;
			}
		}
	}

/** READ EIGENVAL **/
        Eband = MakeArry(2, p.nkpts, p.nbands);
        Eband = ReadEIGENVAL(p, n);

/** WRITE eigen.bd **/
	WriteBanddat(p, kptsco, Eband);

/** WRITE Eigenbp.in & tic.dat **/
	WriteGNUscript(p, hsptsco);

	free(dkpts);
	free(hsptsco);
	free(kptsco);
	free(Eband);

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

double*** ReadEIGENVAL (PARA p, NZW n) {

	FILE *fp;

	int i, j;
	char buffer[512];

	double ***out;
	double tmp;


	if((fp = fopen("EIGENVAL", "r")) == NULL) {
		printf("EIGENVAL is not found.\n");
		return 0;
	}

	out = MakeArry(2, p.nkpts, p.nbands);
	for(i=0; i<6; i++) fgets(buffer, 512, fp);
	if(n.hf == 1) {
		for(i=0; i<n.nnzw; i++) {
			for(j=0; j<p.nbands+2; j++) fgets(buffer, 512, fp);
		}
	}
	if(p.soc==0 && p.ispin==2) {
		for(i=0; i<p.nkpts; i++) {
			for(j=0; j<2; j++) fgets(buffer, 512, fp);
			for(j=0; j<p.nbands; j++) fscanf(fp, "%s %lf %lf %lf %lf", buffer, &out[0][i][j], &out[1][i][j], &tmp, &tmp);
			fgets(buffer, 512, fp);
		}
	}
	else if(p.soc==1 || p.ispin==1) {
		for(i=0; i<p.nkpts; i++) {
			for(j=0; j<2; j++) fgets(buffer, 512, fp);
			for(j=0; j<p.nbands; j++) fscanf(fp, "%s %lf %lf", buffer, &out[0][i][j], &tmp);
			fgets(buffer, 512, fp);
		}
	}

	fclose(fp);

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

double* CalcDKPTS (int nlines, double ***kpts) {

	int i;
	double *out;


	out = (double*)malloc(sizeof(double)*(nlines+1));

        for(i=0; i<nlines; i++) out[i] = sqrt((pow((kpts[i][0][0]-kpts[i][1][0]), 2.0) + pow((kpts[i][0][1]-kpts[i][1][1]), 2.0) + pow((kpts[i][0][2]-kpts[i][1][2]), 2.0)));
	
	return out;

}

void WriteBanddat (PARA p, double *kptsco, double ***Eband) {
	
	FILE *fp;

	int i, j;


	fp = fopen("eigen.bd", "w");

	if(p.soc==0 && p.ispin==2) {
		for(i=0; i<p.nbands; i++) {
			for(j=0; j<p.nkpts; j++) {
				fprintf(fp, "%14.10f %14.10f %14.10f\n", kptsco[j], Eband[0][j][i]-p.Ef, Eband[1][j][i]-p.Ef);
			}
			fprintf(fp, "\n");
		}
	}
	else if(p.soc==1 || p.ispin==1) {
		for(i=0; i<p.nbands; i++) {
			for(j=0; j<p.nkpts; j++) {
				fprintf(fp, "%14.10f %14.10f\n", kptsco[j], Eband[0][j][i]-p.Ef);
			}
			fprintf(fp, "\n");
		}
	}

	fclose(fp);

}

void WriteGNUscript (PARA p, double *hsptsco) {

	FILE *eigen2bp, *tic;

	int i, j;


	if(p.nlines != 1) tic = fopen("tic.dat", "w");
	eigen2bp = fopen("Eigenbp.in", "w");

	fprintf(eigen2bp, "set size 0.5,0.5\n");
	fprintf(eigen2bp, "set terminal postscript enhanced 12\n");
	fprintf(eigen2bp, "set output 'band.eps'\n");
	fprintf(eigen2bp, "set nokey\n");
	fprintf(eigen2bp, "set ylabel 'Energy (eV)'\n");
	fprintf(eigen2bp, "set xrange [0:%f]\n", hsptsco[p.nlines]);
	fprintf(eigen2bp, "set yrange [%f:%f]\n", p.Emin, p.Emax);
	fprintf(eigen2bp, "set ytics %f,%f,%f\n", p.Emin, (p.Emax-p.Emin)/2, p.Emax);
	fprintf(eigen2bp, "set pointsize 0.2\n");
	fprintf(eigen2bp, "set bmargin 3\n");
	fprintf(eigen2bp, "set xzeroaxis\n");
	for(i=0; i<p.nlines+1; i++) {
		if(p.hspts[i]=='G' && i==0)
			fprintf(eigen2bp, "set xtics (""\"{/Symbol %c}\""" %f,", p.hspts[i], hsptsco[i]);
		if(p.hspts[i]!='G' && i==0)
			fprintf(eigen2bp, "set xtics (""\"%c\""" %f,", p.hspts[i], hsptsco[i]);
		if(p.hspts[i]=='G' && i!=0 && i!=p.nlines)
			fprintf(eigen2bp, " ""\"{/Symbol %c}\""" %f,", p.hspts[i], hsptsco[i]);
		if(p.hspts[i]!='G' && i!=0 && i!=p.nlines)
			fprintf(eigen2bp, " ""\"%c\""" %f,", p.hspts[i], hsptsco[i]);
		if(p.hspts[i]=='G' && i==p.nlines)
			fprintf(eigen2bp, " ""\"{/Symbol %c}\""" %f)\n", p.hspts[i], hsptsco[i]);
		if(p.hspts[i]!='G' && i==p.nlines)
			fprintf(eigen2bp, " ""\"%c\""" %f)\n", p.hspts[i], hsptsco[i]);
	}
	if(p.soc==0 && p.ispin==2) {
		fprintf(eigen2bp, "plot	'eigen.bd' using 1:2 with lines lc rgb ""\"black\""" ,\\\n");
		if(p.nlines!=1) fprintf(eigen2bp, "	'eigen.bd' using 1:3 with lines lc rgb ""\"black\""" ,\\\n");
		else if(p.nlines==1) fprintf(eigen2bp, "	'eigen.bd' using 1:3 with lines lc rgb ""\"black\"""");
	}
	else if(p.soc==1 || p.ispin==1) {
		if(p.nlines!=1) fprintf(eigen2bp, "plot	'eigen.bd' using 1:2 with lines lc rgb ""\"black\""" ,\\\n");
		else if(p.nlines==1) fprintf(eigen2bp, "plot	'eigen.bd' using 1:2 with lines lc rgb ""\"black\"""");
	}
	if(p.nlines != 1) {
		fprintf(eigen2bp, "	'tic.dat' using 1:2 with lines\n");
		for(i=0; i<p.nlines-1; i++) for(j=-50; j<51; j++) fprintf(tic, "%f  %d\n", hsptsco[i+1], j);
	}

	fclose(eigen2bp);
	if(p.nlines != 1) fclose(tic);

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
