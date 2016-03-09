#include "s3_reconmerge.h"
/*
void S3_ReconMerge(hls::stream< t_nufft_output_complex > H[reconC],
					hls::stream< t_nufft_output_complex > L0[reconC],
				hls::stream< t_nufft_output_complex > LA[reconC],
				hls::stream< t_nufft_output_complex > LP[reconC],
		        hls::stream< t_recon_complex>  imDFTOut[reconC])
*/
#define reconC    3

int main() {
	hls::stream< t_nufft_output_complex > H;
	hls::stream< t_nufft_output_complex > L0;
	hls::stream< t_nufft_output_complex > LA;
	hls::stream< t_nufft_output_complex > LP;
	hls::stream< t_recon_complex>  imDFTOut;
	for(int k=0;k<8;k++) {
		for(int i=0;i<512;i++) {
				t_nufft_output_complex v ( rand()%100, rand()%100);
				H.write(v);
				L0.write(v);
			}
		for(int i=0;i<480;i++) {
				t_nufft_output_complex v ( rand()%100, rand()%100);
				LA.write(v);
			}
		for(int i=0;i<16;i++) {
				t_nufft_output_complex v ( rand()%100, rand()%100);
				LP.write(v);
		}
	}
	S3_ReconMerge( H, L0, LA, LP, imDFTOut);
	S3_ReconMerge( H, L0, LA, LP, imDFTOut);
	S3_ReconMerge( H, L0, LA, LP, imDFTOut);
	S3_ReconMerge( H, L0, LA, LP, imDFTOut);
	FILE * fo = fopen("pyr_out.txt", " wb");
	for(int i=0;i<512*2*4;i++)
	//for(int i=0;i<2992;i++)
	//while(!pyrFilOut.empty())
	{
		t_recon_complex val = imDFTOut.read();
		//std::cout << val << std::endl;
		fprintf(fo, "%.8f %.8f\n", val.real().to_float(), val.imag().to_float());
		//printf( "%.8f  %.8f\n", val.real().to_float(), val.imag().to_float());
		std::cout << val <<  std::endl;
	}
	fclose(fo);

return 0;
}
