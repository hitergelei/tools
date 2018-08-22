#include <stdio.h>
#include <stdlib.h>

#define bl 1.0
#define cos30 0.866025404

int main(int argc, char* argv[]) {

	FILE *poscar;

	int i;
	double tmp;

	int n, etype;
	double a[3], hul, hw;


        if (argc != 5) {
                printf("Usage: nrgeo <length> <vacuum> <# of hexagon rings(size)> <zigzag(0)/armchair(1)>\n");
		printf("e.g.) nrgeo 1.42 15 1.5 0\n");
                return 0;
        }
	else {
		a[0] = atof(argv[1]);
//		a[1] = atof(argv[2]);
		a[2] = atof(argv[2]);
		tmp = atof(argv[3]);
		etype = atoi(argv[4]);
	}
	if(etype == 0) {
		n = (int)tmp;
		hul = (2*bl+(3*n-2)*(a[0]/2)+a[2]);
		hw = (2*bl+(3*n-2)*(a[0]/2));
	}
	else if(etype == 1) {
		n = (int)(2*tmp);
		hul = ((2*bl+n*a[0])*cos30+a[2]);
		hw = ((2*bl+n*a[0])*cos30);
	}
	

	poscar = fopen("POSCAR", "w");

	fprintf(poscar, "nanoribbon\n");
	fprintf(poscar, "  1.0\n");
	fprintf(poscar, "   %16.12f %16.12f %16.12f\n", (double)((int)(hul+0.9)), 0.0, 0.0);
	if(etype == 0) {
		fprintf(poscar, "   %16.12f %16.12f %16.12f\n", 0.0, 2*a[0]*cos30, 0.0);
	}
	if(etype == 1) {
		fprintf(poscar, "   %16.12f %16.12f %16.12f\n", 0.0, 3*a[0], 0.0);
	}
	fprintf(poscar, "   %16.12f %16.12f %16.12f\n", 0.0, 0.0, a[2]);
	fprintf(poscar, "   C   H\n");
	if(etype == 0 ) {
		fprintf(poscar, "%4d%4d\n", 2*n, 2);
	}
	else if(etype == 1) {
		fprintf(poscar, "%4d%4d\n", 2*(n+1), 4);
	}
	fprintf(poscar, "Cartesian\n");
	hul = ((double)((int)(hul+0.9)) - hw)/2.0;
	if(etype == 0) {
		for(i=0;i<(n+1)/2;i++) {
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+bl+3*i*a[0], a[0]*cos30*3/2, a[2]/2);
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+bl+a[0]/2+3*i*a[0], a[0]*cos30*1/2, a[2]/2);
		}
		for(i=0;i<n/2;i++) {
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+bl+a[0]*3/2+3*i*a[0], a[0]*cos30*1/2, a[2]/2);
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+bl+a[0]*2+3*i*a[0], a[0]*cos30*3/2, a[2]/2);
		}
		fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul, a[0]*cos30*3/2, a[2]/2);
		if(n%2 == 0) {
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+2*bl+a[0]*(3*n-2)/2, a[0]*cos30*3/2, a[2]/2);
		}
		else if(n%2 == 1) {
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+2*bl+a[0]*(3*n-2)/2, a[0]*cos30*1/2, a[2]/2);
		}
	}
	else if(etype == 1) {
		if(n%2 == 0) {
			for(i=0;i<(n/2);i++) {
				fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(bl+2*i*a[0])*cos30, a[0], a[2]/2);
				fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(bl+2*i*a[0])*cos30, 2*a[0], a[2]/2);
				fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(bl+2*i*a[0]+a[0])*cos30, a[0]/2, a[2]/2);
				fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(bl+2*i*a[0]+a[0])*cos30, 5*a[0]/2, a[2]/2);
			}
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(bl+2*i*a[0])*cos30, a[0], a[2]/2);
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(bl+2*i*a[0])*cos30, 2*a[0], a[2]/2);
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul, a[0]-bl/2, a[2]/2);
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul, 2*a[0]+bl/2, a[2]/2);
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(2*bl+2*i*a[0])*cos30, a[0]-bl/2, a[2]/2);
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(2*bl+2*i*a[0])*cos30, 2*a[0]+bl/2, a[2]/2);
		}
		else if(n%2 == 1) {
			for(i=0;i<(n/2)+1;i++) {
				fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(bl+2*i*a[0])*cos30, a[0], a[2]/2);
				fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(bl+2*i*a[0])*cos30, 2*a[0], a[2]/2);
				fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(bl+2*i*a[0]+a[0])*cos30, a[0]/2, a[2]/2);
				fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(bl+2*i*a[0]+a[0])*cos30, 5*a[0]/2, a[2]/2);
			}
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul, a[0]-bl/2, a[2]/2);
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul, 2*a[0]+bl/2, a[2]/2);
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(2*bl+2*i*a[0]-a[0])*cos30, a[0]/2+bl/2, a[2]/2);
			fprintf(poscar, "   %16.12f %16.12f %16.12f\n", hul+(2*bl+2*i*a[0]-a[0])*cos30, 5*a[0]/2-bl/2, a[2]/2);
		}
	}

	fclose(poscar);

	return 0;

}
