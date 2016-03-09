
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
	t_input_complex sigRef[NLEN];
	t_disp_scalar prealign[NLEN];
	t_output_complex cmp[NLEN];
	FILE * fo = fopen("signalLeftPic.txt", " wb");
	FILE * fo2 = fopen("Disparity.txt", " wb");
	const int Kset = 7;
	const int limits[] = { 512, 512,256,128,64,32,16};
	const int climits[] = { 0, 512, 1024,1280,1408,1472,1504};
	const int llimits[] = { 9,9,8,7,6,5,4};
		int l = 0;
		int i = 0;
	  for(int x=0;x<NLEN;x++ ){
		  int nlimit = limits[l];
		  int half =nlimit>>1;
		  int idx;
		  int limit = limits[l];
		  int cidx  = climits[l];
		  	  if(i<half) idx = i;
		  	  else idx = 512 - nlimit + i;
		  	sig[x].real()    = pyrLbasicI[idx * 2 + 0];
		  	sig[x].imag() 	 = pyrLbasicI[idx * 2 + 1];
		  	sigRef[x].real()    = pyrRdetailI[idx * 2 + 0];
		  	sigRef[x].imag() = pyrRdetailI[idx * 2 + 1];
		  	fprintf(fo, "%.8f %.8f +Disp  %.8f \n", sig[i].real().to_float(), sig[i].imag().to_float(),pyrLbasicI[i * 2 + 0]);
		  	if (x<512)prealign[x]  = prealignI[x];
		  	fprintf(fo2, "%.8f  \n", prealign[i].to_float());
		  		 }

	fclose(fo);
	fclose(fo2);

	fo = fopen("output.txt", " wb");
	int nL = 512;
	int nLExp = 1024;
	int width = 1024;
	cmpy_complex_top(sig, sigRef, prealign, cmp, nL, nLExp, nL, -1.0*nL/width);
	for(int i=0;i<1520;i++ ){
		fprintf(fo,"%.8f  %.8f \n ", cmp[i].real().to_float(), cmp[i].imag().to_float(), cmpRI[i], cmpII[i]);
		printf("%.8f %.8f %.8f %.8f\n ", cmp[i].real().to_float(), cmp[i].imag().to_float(), cmpRI[i], cmpII[i]);
	}
	fclose(fo);
	return 0;
}
