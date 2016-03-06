
#include "cmpy_complex.h"
#include <stdio.h>
int main() {

	float pyrLbasicI[8192*2];
	float pyrRdetailI[8192*2];
	float prealignI[8192];

	float cmpRI[8192];
	float cmpII[8192];

	FILE *fi = NULL;


	fi = fopen("pyrLbasic.dat", "rb");
	fread(pyrLbasicI, sizeof(float)*2, 8192, fi);
	fclose(fi);

	fi = fopen("pyrRdetail.dat", "rb");
	fread(pyrRdetailI, sizeof(float)*2, 8192, fi);
	fclose(fi);

	fi = fopen("prealign.dat", "rb");
	fread(prealignI, sizeof(float), 8192, fi);
	fclose(fi);

	fi = fopen("pyrCompIm.dat", "rb");
	fread(cmpII, sizeof(float), 8192, fi);
	fclose(fi);

	fi = fopen("pyrCompRe.dat", "rb");
	fread(cmpRI, sizeof(float), 8192, fi);
	fclose(fi);

	printf("Read data %x\n", fi);
	t_input_complex sig[NLEN];
	t_input_complex sigRef[NLEN * 2];
	t_disp_scalar prealign[NLEN];
	t_output_complex cmp[NLEN];
	FILE * fo = fopen("signalLeftPic.txt", " wb");
	FILE * fo2 = fopen("Disparity.txt", " wb");
	for(int i=0;i<NLEN;i++) {
		sig[i].real()    = pyrLbasicI[i * 2 + 0];
		sig[i].imag() = pyrLbasicI[i * 2 + 1];
		fprintf(fo, "%.8f %.8f +Disp  %.8f \n", sig[i].real().to_float(), sig[i].imag().to_float(),pyrLbasicI[i * 2 + 0]);
		//sig[i].last = 0;
		prealign[i]  = prealignI[i];
		fprintf(fo2, "%.8f  \n", prealign[i].to_float());
	}
	fclose(fo);
	fclose(fo2);
	fo = fopen("signalRightPic.txt", " wb");
	for(int i=0;i<NLEN*2;i++) {
		sigRef[i].real()    = pyrRdetailI[i * 2 + 0];
		sigRef[i].imag() = pyrRdetailI[i * 2 + 1];
		fprintf(fo, "%.8f %.8f\n", sigRef[i].real().to_float(), sigRef[i].imag().to_float());
	}
	fclose(fo);
	fo = fopen("output.txt", " wb");
	int nL = 512;
	int nLExp = 1024;
	int width = 1024;
	cmpy_complex_top(sig, sigRef, prealign, cmp, nL, nLExp, nL, -1.0*nL/width);
	for(int i=0;i<512;i++ ){
		fprintf(fo,"%.8f  %.8f \n ", cmp[i].real().to_float(), cmp[i].imag().to_float(), cmpRI[i], cmpII[i]);
		printf("%.8f %.8f %.8f %.8f\n ", cmp[i].real().to_float(), cmp[i].imag().to_float(), cmpRI[i], cmpII[i]);
	}
	fclose(fo);
	return 0;
}
