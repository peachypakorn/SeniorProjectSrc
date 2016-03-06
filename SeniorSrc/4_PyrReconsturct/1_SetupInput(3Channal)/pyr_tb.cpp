#include "s0_nufftinput.h"

typedef    unsigned short dataInF[7];
int main() {
	dataInF pyr[1520*2];
	FILE *fi = fopen("pyrIn112b.dat", "r");
	fread(pyr, 1520*2, sizeof(char) * 14, fi);
	fclose(fi);



//printf("%d", sizeof(packed));
	hls::stream< packed >    src_axi;
	// Writing data
	for(int i=0;i<1520;i++) {
		packed val;
		for(int v=0;v<7;v++)
			val.range((v+1)*16 -1, v*16) = pyr[i][v];
		src_axi.write( val);
	}

	hls::stream<t_input_complex> S0_pyr[4][3];
	hls::stream<t_disp_scalar>   S0_disp[4];
	nufft_input_interface(src_axi,S0_pyr,S0_disp);

	FILE * fo = fopen("OutputOfTestInput.txt", " wb");
	fprintf(fo,"Data Real Disp(Phase");
		for(int i=0;i<1520;i++)
		{
			t_input_complex val = S0_pyr[0][0].read();
			 S0_pyr[0][1].read();
			 S0_pyr[0][2].read();
			 S0_pyr[1][0].read();
			 S0_pyr[1][1].read();
			 S0_pyr[1][2].read();
			 S0_pyr[2][0].read();
			 S0_pyr[2][1].read();
			 			 S0_pyr[2][2].read();
			 			 S0_pyr[3][0].read();
			 			 S0_pyr[3][1].read();
			 			 S0_pyr[3][2].read();

			fprintf(fo, "%.8f %.8f  Disp:   ", val.real().to_float(), val.imag().to_float());
			t_disp_scalar d1 ;
			for (int i = 0; i < 4; ++i) {
				 d1 = S0_disp[i].read();
			fprintf(fo, "%.8f ", d1.to_float());
			}
			fprintf(fo,"\n");

		}
		fclose(fo);
	printf("Done!!!!!!!!!!!!!!!!!!!!!!!!\n");
		return 0;
}
