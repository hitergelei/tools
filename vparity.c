#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define	NG		200
#define	TPLWV	10000
#define	PI		3.14159265358979323846264338327950
#define	ATOMICA		0.529177249
#define RYTOEV		13.605826
#define TINY		(1.0e-10)

// check_parity
int	check_parity(double cpwv[][2], int npwv, int gvec[][3], double kpt[3]);
// matrix, vector operators
double	dot(double v1[3], double v2[3]);
double	length(double v[3]);
int	cross(double v1[3], double v2[3], double r[3]);
double	det(double v[3][3]);
// generate reciprocal lattice vectors
double	gen_bvec(double avec[3][3], double bvec[3][3]);
// calculate the length
double	len_k(double c[3], double bvec[3][3]);
// read fermi level from outcar
double read_outcar_fermi(void);

typedef struct {
	double	nkpt;
	double	nband;
	double	enmax;
	double	a[3][3];
} WDES;

typedef struct {
	double	npl;
	double	vec[3];
	double	eig[][3];
} KPT;

typedef struct {
	int	nkpt;
	int	nband;
	double	kpt[8][3];
	char	v[4];
	char	pi[8];
	char	parity[][8];
} PRT;

int	main(void)
{
	long	rec, ispin, tag;
	double	*pbuff;
	double	buff[18000];
	float	buff_[9000];
	WDES	data;
	KPT	**pkpt;
	long	lsize;
	FILE*	fwave;
	size_t	readb;
	int	i, j, k, n, p;
	//
	double	ecut;
	double	coe[3];
	double	vol;
	double	b[3][3];
	int	mindx[NG];
	int	npwv, nx, ny, nz;
	int	gvec[TPLWV][3];
	double	cpwv[TPLWV][2];
	//
	int	parity, ptot =1;
	double	fermi;

	fermi = read_outcar_fermi();

	fwave = fopen("WAVECAR", "rb");
	if (fwave == NULL)
	{
		printf("Failed to open WAVECAR\n");
		return 0;
	}

	fseek( fwave, 0, SEEK_END);
	lsize = ftell( fwave);
	printf("File size %ld \n",lsize);
	// read initial data (rec size, ispin, tag)
	fseek( fwave, 0, SEEK_SET);
	readb = fread(buff, sizeof(double), 3, fwave);

	rec = (long)buff[0];
	ispin = (long)buff[1];
	tag = (long)buff[2];
	printf("rec = %ld, ispin = %ld, tag = %ld\n",
		rec, ispin, tag);
	// read nkpt, nband, emax
	fseek( fwave, rec, SEEK_SET);
	readb = fread(&data, sizeof(data), 1, fwave);

	ecut = data.enmax/RYTOEV / ATOMICA / ATOMICA;
	vol = gen_bvec(data.a, b);
	printf("nkpt = %f, nband = %f, emax = %f\n",
		data.nkpt, data.nband, data.enmax);
	printf("Primitive vector (vol = %f): \n", vol);
	printf("%f %f %f\n",data.a[0][0], data.a[0][1], data.a[0][2]);
	printf("%f %f %f\n",data.a[1][0], data.a[1][1], data.a[1][2]);
	printf("%f %f %f\n",data.a[2][0], data.a[2][1], data.a[2][2]);
	printf("Reciprocal vector : \n");
	printf("%f %f %f\n",b[0][0], b[0][1], b[0][2]);
	printf("%f %f %f\n",b[1][0], b[1][1], b[1][2]);
	printf("%f %f %f\n\n",b[2][0], b[2][1], b[2][2]);

	//pbuff = (double*)calloc(rec, sizeof(char));
	pkpt = (KPT**)calloc((int)data.nkpt, sizeof(KPT*));
	lsize = 4+3*(int)data.nband; // sizeof KPT
	for (i=0;i<(int)data.nkpt;i++)
		pkpt[i] = (KPT*)calloc(lsize, sizeof(double));

	// set m index : ( 0, 1, 2, ... , NG/2, -(NG/2 -1), ..., -1 )
	for (i=0;i<=NG/2;i++)
		mindx[i] = i;
	for (i=NG/2+1;i<NG;i++)
		mindx[i] = i-NG;

	n = 1;
	for (i=0;i<(int)ispin;i++)
	for (j=0;j<(int)data.nkpt;j++)
	{
		// read kpoint information
		n++;
		fseek( fwave, rec*n, SEEK_SET);
		readb = fread(pkpt[j], sizeof(double)*lsize, 1, fwave);

		printf("kpt %3d (%e,%e,%e)  npl (%.0f)\n", j+1,
			pkpt[j]->vec[0], pkpt[j]->vec[1], pkpt[j]->vec[2],
			pkpt[j]->npl);
//		for (i=0;i<(int)data.nband;i++)
//			printf("%f %f %f\n", pkpt[j]->eig[i][0],
//				pkpt[j]->eig[i][1],
//				pkpt[j]->eig[i][2]);

		npwv = 0;
		for(nz=0;nz<NG;nz++)
		for(ny=0;ny<NG;ny++)
		for(nx=0;nx<NG;nx++)
		{
			coe[0] = (double)mindx[nx] + pkpt[j]->vec[0];
			coe[1] = (double)mindx[ny] + pkpt[j]->vec[1];
			coe[2] = (double)mindx[nz] + pkpt[j]->vec[2];
			if (len_k(coe, b) < ecut)
			{
				gvec[npwv][0] = mindx[nx];
				gvec[npwv][1] = mindx[ny];
				gvec[npwv][2] = mindx[nz];
				npwv++;
			}
		}
		if (pkpt[j]->npl != npwv)
			printf(" !!!!! npwv(%d) is different from pkpt[j]->npl(%.0f)\n", npwv, pkpt[j]->npl);

		ptot = 1;
		//for (k=0;k<2*(int)data.nband;k++)
		for (k=0;k<(int)data.nband;k++)
		{
			// read wavefunction info for each band
			n++;
			fseek( fwave, rec*n, SEEK_SET);
			readb = fread(buff_, sizeof(buff_), 1, fwave);

			for (p=0;p<npwv;p++)
			{
				cpwv[p][0] = (double)buff_[2*p];
				cpwv[p][1] = (double)buff_[2*p+1];
		//	printf("cpwv[%3d] = [% d,% d,% d] (% e,% e)\n", 
		//		p, gvec[p][0], gvec[p][1], gvec[p][2], 
		//		cpwv[p][0], cpwv[p][1]);
			}
			parity = check_parity(cpwv, npwv, gvec, pkpt[j]->vec);
			printf("band %3d : %c   % .5f\n", k,
			(parity==1)? '+' : (parity==-1)? '-' : 'o',
			pkpt[j]->eig[k][0]);
			if (pkpt[j]->eig[k][0] < fermi)
			{
				ptot *= parity;
			if (k!=(int)data.nband-1 &&
				fermi < pkpt[j]->eig[k+1][0])
			printf("---------- parity [%c] ----------\n",
			(ptot==1)? '+' : (ptot==-1)? '-' : 'o');
			}
		}
	}

	for (i=0;i<(int)data.nkpt;i++)
		free(pkpt[i]);
	free(pkpt);

	fclose(fwave);

	return 1;
}

int	check_parity(double cpwv[][2], int npwv, int gvec[][3], double kpt[3])
{
	int	i, j, n[3];
	int	pnew[2], parity = 0;
	double	diff[2], acpwv[2][2];
	n[0] = (kpt[0]==0.0) ? 0 : 1;
	n[1] = (kpt[1]==0.0) ? 0 : 1;
	n[2] = (kpt[2]==0.0) ? 0 : 1;
	for (i=1;i<npwv;i++)
	for (j=i;j<npwv;j++)
	{
		if (gvec[i][0] == -gvec[j][0] -n[0] &&
		    gvec[i][1] == -gvec[j][1] -n[1] &&
		    gvec[i][2] == -gvec[j][2] -n[2])
		{
			pnew[0] = (cpwv[i][0]*cpwv[j][0]>=0) ? 1 : -1;
			pnew[1] = (cpwv[i][1]*cpwv[j][1]>=0) ? 1 : -1;
			acpwv[0][0] = fabs(cpwv[i][0]);
			acpwv[0][1] = fabs(cpwv[i][1]);
			acpwv[1][0] = fabs(cpwv[j][0]);
			acpwv[1][1] = fabs(cpwv[j][1]);

			if (pnew[0]==1) 
				diff[0] = fabs(cpwv[i][0]-cpwv[j][0]);
			else	diff[0] = fabs(cpwv[i][0]+cpwv[j][0]);
			if (pnew[1]==1) 
				diff[1] = fabs(cpwv[i][1]-cpwv[j][1]);
			else	diff[1] = fabs(cpwv[i][1]+cpwv[j][1]);

			if (acpwv[0][0] > 0.00001 &&
			    acpwv[0][1] > 0.00001 &&
			    acpwv[1][1] > 0.00001 &&
			    acpwv[1][1] > 0.00001 &&
			        diff[0] < 0.00001 &&
			        diff[1] < 0.00001)
			{
			if (pnew[0] != pnew[1])
			{
printf("[%4d] (% d,% d,% d) (% e,% e)\n[%4d] (% d,% d,% d) (% e,% e)",
i, gvec[i][0], gvec[i][1], gvec[i][2], cpwv[i][0], cpwv[i][1],
j, gvec[j][0], gvec[j][1], gvec[j][2], cpwv[j][0], cpwv[j][1]);
printf("  Difference (% e,% e)\n", diff[0], diff[1]);
			return 0;
			}
			if (parity == 0) // set parity
				parity = pnew[0];
			else if (parity != pnew[0])
			{
printf("[%4d] (% d,% d,% d) (% e,% e)\n[%4d] (% d,% d,% d) (% e,% e)",
i, gvec[i][0], gvec[i][1], gvec[i][2], cpwv[i][0], cpwv[i][1],
j, gvec[j][0], gvec[j][1], gvec[j][2], cpwv[j][0], cpwv[j][1]);
printf("  Difference (% e,% e)\n", diff[0], diff[1]);
			return 0;
			}
			}
		} 
	}
	return parity;
}

double	dot(double v1[3], double v2[3])
{
	int	i;
	double	ret=0.0;
	for (i=0;i<3;i++)
		ret += v1[i] * v2[i];
	return ret;
}

double	length(double v[3])
{
	return sqrt(dot(v,v));
}

int	cross(double v1[3], double v2[3], double r[3])
{
	int	i;
	for (i=0;i<3;i++)
	{
		r[i] = v1[(i+1)%3] * v2[(i+2)%3]
		     - v1[(i+2)%3] * v2[(i+1)%3];
	}
	return 1;
}

double	det(double v[3][3])
{
	double r[3];
	cross(v[1], v[2], r);
	return dot(v[0], r);
}

double	gen_bvec(double avec[3][3], double bvec[3][3])
{
	int	i, j;
	double	vol;

	vol = det(avec);
	for (i=0;i<3;i++)
		cross(avec[(i+1)%3],avec[(i+2)%3],bvec[i]);
	for (i=0;i<3;i++)
	for (j=0;j<3;j++)
		bvec[i][j]*=(2.0*PI/vol);
		bvec[i][j]/=vol;
	return vol;
}

double	len_k(double c[3], double bvec[3][3])
{
	double	v[3] = {0.0,0.0,0.0};
	int	i,j;

	for (i=0;i<3;i++)
	for (j=0;j<3;j++)
		v[i] += c[j]*bvec[j][i];

	return dot(v,v);
}

double read_outcar_fermi(void)  //read fermi level
{
        FILE	*fp;
        int	i;
	char	buff[512];
        double	fermi;
        if((fp=fopen("OUTCAR","r"))==NULL)fp=fopen("outcar","r");
        fscanf(fp,"%s",buff);
        while((buff[0]!='E')||(buff[1]!='-'))
        {
                fgets(buff,512,fp);
                fscanf(fp,"%s",buff);
        }
        fscanf(fp,"%s %lf",buff,&fermi);
        return fermi;
}
