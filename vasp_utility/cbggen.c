#include <stdio.h>
#include <stdlib.h>

#define bl 1.5
#define cos30 0.866025404

int main(int argc, char* argv[]) {

	FILE *output;

	int i, n, h;
	double a, tl;


	if (argc != 4) {
		printf("Usage: cbggen <bond length> <# of hexagon rings(size)> <0/1>\n");
		printf("e.g.) cbggen 3 1.5 0\n");
		return 0;
        }
	else {
		a = atof(argv[1]);
		n = (int)2*atof(argv[2]) + 1;
		h = atoi(argv[3]);
	}
	tl = a*(n-1)*cos30;

	output = fopen("POSCAR", "w");

	fprintf(output, "CBG-like structure\n");
	fprintf(output, "  1.0\n");
	fprintf(output, "   %16.12f %16.12f %16.12f\n", 1.5*a*(n-1)*cos30, -(n-1)*a*cos30*cos30, 0.0);
	fprintf(output, "   %16.12f %16.12f %16.12f\n", 1.5*a*(n-1)*cos30, (n-1)*a*cos30*cos30, 0.0);
	fprintf(output, "   %16.12f %16.12f %16.12f\n", 0.0, 0.0, 3*a);
	if(h == 0) fprintf(output, "  Bi\n");
	else if(h == 1) fprintf(output, "  Bi   H\n");
	if(h == 0) fprintf(output, "%4d\n", 10+(n-3)*6);
	else if(h == 1) fprintf(output, "%4d%4d\n", 10+(n-3)*6, (n-2)*6);
	fprintf(output, "Cartesian\n");
	for(i=0; i<n; i++) {
		if(i%2 == 0) {
			fprintf(output, "   %16.12f %16.12f %16.12f\n", (i*a*cos30+tl), 0.0, 0.5*a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n", (i*a*cos30+tl), 0.0, 2.5*a);
		}
		else {
			fprintf(output, "   %16.12f %16.12f %16.12f\n", (i*a*cos30+tl), 0.0, a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n", (i*a*cos30+tl), 0.0, 2*a);
		}
	}
	for(i=0; i<n-2; i++) {
		if((i%2 == 0) && (n%2 != 0)) {
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), (i+1)*a*cos30*cos30, a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), (i+1)*a*cos30*cos30, 2*a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), -(i+1)*a*cos30*cos30, a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), -(i+1)*a*cos30*cos30, 2*a);
		}
		else if((i%2 == 0) && (n%2 == 0)) {
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), (i+1)*a*cos30*cos30, 0.5*a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), (i+1)*a*cos30*cos30, 2.5*a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), -(i+1)*a*cos30*cos30, 0.5*a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), -(i+1)*a*cos30*cos30, 2.5*a);
		}
		else if((i%2 != 0) && (n%2 != 0)) {
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), (i+1)*a*cos30*cos30, 0.5*a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), (i+1)*a*cos30*cos30, 2.5*a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), -(i+1)*a*cos30*cos30, 0.5*a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), -(i+1)*a*cos30*cos30, 2.5*a);
		}
		else if((i%2 != 0) && (n%2 == 0)) {
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), (i+1)*a*cos30*cos30, a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), (i+1)*a*cos30*cos30, 2*a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), -(i+1)*a*cos30*cos30, a);
			fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+tl), -(i+1)*a*cos30*cos30, 2*a);
		}
	}

	if(h == 1) {
		for(i=0; i<n; i++) {
			if(i%2 == 0 && i != 0 && i != n-1) {
				fprintf(output, "   %16.12f %16.12f %16.12f\n", (i*a*cos30+tl), -bl, 0.5*a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n", (i*a*cos30+tl), bl, 2.5*a);
			}
			else if(i != 0 && i != n-1) {
				fprintf(output, "   %16.12f %16.12f %16.12f\n", (i*a*cos30+tl), bl, a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n", (i*a*cos30+tl),-bl, 2*a);
			}
		}
		for(i=0; i<n-2; i++) {
			if((i%2 == 0) && (n%2 != 0)) {
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+bl*cos30+tl), ((i+1)*a*cos30*cos30-bl*0.5), a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30-bl*cos30+tl), ((i+1)*a*cos30*cos30+bl*0.5), 2*a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30-bl*cos30+tl), -((i+1)*a*cos30*cos30+bl*0.5), a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+bl*cos30+tl), -((i+1)*a*cos30*cos30-bl*0.5), 2*a);
			}
			else if((i%2 == 0) && (n%2 == 0)) {
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30-bl*cos30+tl), ((i+1)*a*cos30*cos30+bl*0.5), 0.5*a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+bl*cos30+tl), ((i+1)*a*cos30*cos30-bl*0.5), 2.5*a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+bl*cos30+tl), -((i+1)*a*cos30*cos30-bl*0.5), 0.5*a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30-bl*cos30+tl), -((i+1)*a*cos30*cos30+bl*0.5), 2.5*a);
			}
			else if((i%2 != 0) && (n%2 != 0)) {
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30-bl*cos30+tl), ((i+1)*a*cos30*cos30+bl*0.5), 0.5*a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+bl*cos30+tl), ((i+1)*a*cos30*cos30-bl*0.5), 2.5*a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+bl*cos30+tl), -((i+1)*a*cos30*cos30-bl*0.5), 0.5*a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30-bl*cos30+tl), -((i+1)*a*cos30*cos30+bl*0.5), 2.5*a);
			}
			else if((i%2 != 0) && (n%2 == 0)) {
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+bl*cos30+tl), ((i+1)*a*cos30*cos30-bl*0.5), a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30-bl*cos30+tl), ((i+1)*a*cos30*cos30+bl*0.5), 2*a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30-bl*cos30+tl), -((i+1)*a*cos30*cos30+bl*0.5), a);
				fprintf(output, "   %16.12f %16.12f %16.12f\n",((n-0.5+i*0.5)*a*cos30+bl*cos30+tl), -((i+1)*a*cos30*cos30-bl*0.5), 2*a);
			}
		}
	}

	fclose(output);

	return 0;

}
