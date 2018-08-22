#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {

	FILE *eigenval;
	FILE *hfeigenval;

	int i;
	char buffer[512];

	int nbands, nnzw;


	if(argc != 2) {
		printf("Usage: eigenHF <# of non-zero weight>\n");
		return 0;
	}
	else {
		if((eigenval = fopen("EIGENVAL", "r")) == NULL) {
			printf("OUTCAR is not found.\n");
			return 0;
		}
		hfeigenval = fopen("EIGENVAL_hf","w");
		nnzw = atoi(argv[1]);
	}

        for(i=0; i<5; i++) {
                fgets(buffer, 512, eigenval);
                fputs(buffer, hfeigenval);
        } // 1st~5th lines

	fscanf(eigenval, "%d", &nbands);
	fprintf(hfeigenval, "   %d", nbands);
        fscanf(eigenval, "%d", &nbands);
        fprintf(hfeigenval, "   %d", nbands-nnzw);
        fscanf(eigenval, "%d", &nbands);
        fprintf(hfeigenval, "   %d", nbands);
	// 6th line

        for(i=0; i<2; i++) fgets(buffer, 512, eigenval);
        for(i=0; i<2; i++) fprintf(hfeigenval, "\n"); // 7th~8th lines
	
	for(i=0; i<(nbands+2)*nnzw; i++) fgets(buffer, 512, eigenval); // non-zero weight

	while(1) {
		fgets(buffer, 512, eigenval);
		if(feof(eigenval)) break;
		fputs(buffer, hfeigenval);
	}

	fclose(eigenval);
	fclose(hfeigenval);

	return 0;

}
