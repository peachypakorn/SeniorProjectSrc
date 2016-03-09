#include "s4_finalrecon.h"
/*
void S3_ReconMerge(hls::stream< t_nufft_output_complex > H[reconC],
					hls::stream< t_nufft_output_complex > L0[reconC],
				hls::stream< t_nufft_output_complex > LA[reconC],
				hls::stream< t_nufft_output_complex > LP[reconC],
		        hls::stream< t_recon_complex>  imDFTOut[reconC])
*/
#define reconC    3

int main() {
	hls::stream< t_recon_complex > DFTin;

	hls::stream< PIXEL_RAW>  sigOut;
	for(int k=0;k<8;k++) {
		for(int i=0;i<512*4;i++) {
			t_recon_complex v ( rand()%100, rand()%100);
				DFTin.write(v);

			}

	}
S4_finalrecon(DFTin,sigOut);
S4_finalrecon(DFTin,sigOut);
S4_finalrecon(DFTin,sigOut);
S4_finalrecon(DFTin,sigOut);

	FILE * fo = fopen("pyr_out.txt", " wb");
	for(int i=0;i<512*4;i++)
	//for(int i=0;i<2992;i++)
	//while(!pyrFilOut.empty())
	{
		PIXEL_RAW val = sigOut.read();
		//std::cout << val << std::endl;

		fprintf(fo, "%d\n",val.to_uint()&255 );
		//printf( "%.8f  %.8f\n", val.real().to_float(), val.imag().to_float());
		std::cout << (val.to_uint()&255) <<  std::endl;
	}
	fclose(fo);

return 0;
}
