#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[]) {

	FILE *outcar, *procar;
	FILE *hfprocar;

	int i, j, k, count;
	char buffer[512];

	int nnzw, n[3];
	char soc;


	if(argc != 2) {
		printf("Usage: procarHF <# of non-zero weight>\n");
		return 0;
	}
	else {
		if((procar=fopen("PROCAR", "r")) == NULL) {
			printf("PROCAR is not found.\n");
			return 0;
		}
		if((outcar=fopen("OUTCAR", "r")) == NULL) {
			printf("OUTCAR is not found.\n");
			return 0;
		}
		hfprocar = fopen("PROCAR_HF", "w");
		nnzw = atoi(argv[1]);
	}

	fscanf(outcar, "%s", buffer);
        while(memcmp(buffer, "LNONCOLLINEAR", strlen("LNONCOLLINEAR")) != 0) fscanf(outcar, "%s", buffer);
	for(i=0; i<2; i++) fscanf(outcar, "%s", &soc);

	fclose(outcar);


	fgets(buffer, 512, procar);
	fputs(buffer, hfprocar); // 1st line

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			fscanf(procar, "%s", buffer);
			fprintf(hfprocar, "%s ", buffer);
		}
		fscanf(procar, "%d", &n[i]);
		fprintf(hfprocar, "   %d    ", n[i]);
	}
        fgets(buffer, 512, procar);
        fprintf(hfprocar, "\n"); // 2nd line
	
	fgets(buffer, 512, procar);
	fputs(buffer, hfprocar); // 3rd line
		
	count = 3;
	if(soc == 'F') {
		for(i=0; i<((n[2]+5)*n[1]+3)*nnzw; i++) {
			fgets(buffer, 512, procar);
			count++;
		}
		while(1) {
			fgets(buffer, 512, procar);
			count++;
			fputs(buffer, hfprocar);
			if(count == ((n[2]+5)*n[1]+3)*n[0]+1) break;
		} // spin up
	
		for(i=0; i<3; i++) {
			fgets(buffer, 512, procar);
		        fputs(buffer, hfprocar);
		}
	        for(i=0; i<((n[2]+5)*n[1]+3)*nnzw;i++) {
	                fgets(buffer, 512, procar);
	                count++;
	        }
	        while(1) {
	                fgets(buffer, 512, procar);
	                count++;
	                fputs(buffer, hfprocar);
	                if(feof(procar)) break;
	        } // spin down
	}

	if(soc == 'T') {
                for(i=0; i<((4*n[2]+8)*n[1]+3)*nnzw; i++) {
                        fgets(buffer, 512, procar);
                        count++;
                }
                while(1) {
                        fgets(buffer, 512, procar);
                        count++;
                        fputs(buffer, hfprocar);
                        if(feof(procar)) break;
                }
	}

	fclose(procar);
	fclose(hfprocar);

	return 0;

}
