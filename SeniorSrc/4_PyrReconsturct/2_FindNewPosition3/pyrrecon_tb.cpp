
#include <math.h>

//#include "nufft_cosfilter.h"
#if 0

#include "pyrrecon.h"

void pyrrecon_tb() {
	float pyr[1520*2];

	FILE * fi = fopen ( "pyr_ifft.bin" , "rb");
	fread(pyr, 1520*2, sizeof(float), fi);
	fclose(fi);

	hls::stream< cmpxDataOut>  pyrFFT[reconC];
	hls::stream< t_recon_complex>  imDFTOut[reconC];
	hls::stream< t_output_scalar>  sigOut[reconC];

//	void pyr_recon(hls::stream< t_input_scalar>  pyrFFT[C],
//			       hls::stream< t_input_scalar>  imDFTOut[C]);

//	void pyr_recon_ifft(hls::stream< t_input_scalar>  DFTOut[C],  hls::stream< t_output_scalar>  sigOut[C]);
	for(int i=0;i<1520;i++) {
		cmpxDataOut val;
		val.real() = ((data_out_t) pyr[i*2]) >> 2 ;
		val.imag() = ((data_out_t) pyr[i*2+1]) >> 2 ;
		//printf("%f %f\n", pyr[i*2], pyr[i*2+1]);
		printf("%f + %fi, ", val.real().to_float(), val.imag().to_float());
		for(int c=0;c<reconC;c++)
			pyrFFT[c].write(val);
	}
	printf("\n");

	//pyr_recon( pyrFFT, imDFTOut);
	//pyr_recon_ifft( imDFTOut, sigOut);
	pyr_recon_combine(pyrFFT, sigOut);

	printf("\n R = [");

	for(int i=0;i<512;i++){
//		cmpxDataOut val = imDFTOut[0].read();
//		printf(":%f %f\n", val.real().to_float(), val.imag().to_float());
		for(int c=0;c<reconC;c++) {
			t_output_scalar val = sigOut[c].read();
			printf("%f, ", val.to_float());
		}
	}
	printf("];\n");


	// Rerun
	for(int i=0;i<1520;i++) {
		cmpxDataOut val;
		val.real() = ((data_out_t) pyr[i*2]) >> 2 ;
		val.imag() = ((data_out_t) pyr[i*2+1]) >> 2 ;
		//printf("%f %f\n", pyr[i*2], pyr[i*2+1]);
		printf("%f + %fi, ", val.real().to_float(), val.imag().to_float());
		for(int c=0;c<reconC;c++)
			pyrFFT[c].write(val);
	}
	printf("\n");

	//pyr_recon( pyrFFT, imDFTOut);
	//pyr_recon_ifft( imDFTOut, sigOut);
	pyr_recon_combine(pyrFFT, sigOut);

	printf("\n R = [");

	for(int i=0;i<512;i++){
//		cmpxDataOut val = imDFTOut[0].read();
//		printf(":%f %f\n", val.real().to_float(), val.imag().to_float());
		for(int c=0;c<reconC;c++) {
			t_output_scalar val = sigOut[c].read();
			printf("%f, ", val.to_float());
		}
	}
	printf("];\n");


}

#include "nufft.h"

// Forward declaration
#define C   3
void nufft_top_t(hls::stream<t_input_complex>  sig[C],
			   hls::stream<t_disp_scalar>    &dispFilter,
			   hls::stream<t_nufft_output_complex> sigStreamOut[C]) ;


void nufft_tb() {

	// DEPRECIATE
	FILE * fi = fopen ( "pyr_out.bin" , "rb");
	printf("%x\n", fi);
	float pyr[2992*2];

	t_input_complex  pyrFFT[512];
	t_disp_scalar  	 disp[512];


	fseek(fi, 2992*2 * 39 * sizeof(float), SEEK_SET);
	fread(pyr, 2929*2, sizeof(float), fi);
	fclose(fi);

	int sidx [] =  {0, 512, 1024+512,  1792+512,  2432+256, 2752+128, 2912+64};
	int ssidx [] = {512, 512,256,128,64,32,16};


	FILE * fo = NULL;
	int l=0;
	fo = fopen("pyr_level_out.txt", "wb");
	for(int i=0;i<ssidx[l];i++) {

		//fprintf(fo, "%.8f %.8f\n", pyr[ (sidx[l] + i)*2], pyr[ (sidx[l] + i)*2 + 1]);

		pyrFFT[i] = t_input_complex(pyr[ (sidx[l] + i)*2], pyr[ (sidx[l] + i)*2 + 1]);
		float  dispF =sin(((float)(i+1))/ssidx[l]*4)*4;
		disp[i] = dispF;
//		printf("%f %f\n", dispF, disp[i].to_float());
		fprintf(fo, "%.8f %.8f %.8f\n", pyrFFT[i].real().to_float(), pyrFFT[i].imag().to_float(), disp[i].to_float());

	}
	fclose(fo);


	hls::stream<t_input_complex>  sig[C];
	hls::stream<t_disp_scalar>    dispFilter;
	hls::stream<t_nufft_output_complex> sigStreamOut[C];

	for(int i=0;i<ssidx[l];i++) {
		for(int c=0;c<C;c++) {
			sig[c].write(pyrFFT[i]);
		}
		dispFilter.write(disp[i]);
	}
	nufft_top_t( sig,
			   dispFilter,
			   sigStreamOut);

	fo = fopen("pyr_disp_out.txt", "wb");
	for(int i=0;i<ssidx[l]*2;i++) {
		t_nufft_output_complex vals[C];
		for(int c=0;c<C;c++) {
			vals[c] = sigStreamOut[c].read();
		}
		fprintf(fo, "%.8f %.8f\n", vals[0].real().to_float(), vals[0].imag().to_float());
	}
	fclose(fo);
}



void nufft_tb3() {
//	FILE * fi = fopen ( "pyr_out.bin" , "rb");
//	printf("%x\n", fi);
//	float pyr[2992*2];
	float pyr[1520*2];
	float dispF[1520];

	FILE * fi = fopen("NUFFTWARP_PYRIN.dat", "rb");
	fseek(fi, 1520*2 * 59 * sizeof(float), SEEK_SET);
	fread(pyr, 1520*2, sizeof(float), fi);
	fclose(fi);

	fi = fopen("NUFFTWARP_DISIN.dat", "rb");
	fseek(fi, 1520 * 59 * sizeof(float), SEEK_SET);
	fread(dispF, 1520, sizeof(float), fi);
	fclose(fi);


	t_input_complex  pyrFFT[512];
	t_disp_scalar  	 disp[512];


	int sidx [] =  {0, 512, 1024+512,  1024+256+512,  1024+256+512+128, 1024+256+512 + 64, 1024+256+512 + 64+32};
	int ssidx [] = {512, 512,256,128,64,32,16};


	FILE * fo = NULL;
	int l=1;
	fo = fopen("pyr_level_out.txt", "wb");
	for(int i=0;i<ssidx[l];i++) {

		//fprintf(fo, "%.8f %.8f\n", pyr[ (sidx[l] + i)*2], pyr[ (sidx[l] + i)*2 + 1]);
		int idx = (sidx[l] + i);

		pyrFFT[i] = t_input_complex(pyr[ idx*2], pyr[ idx*2 + 1]);
//		float  dispF =sin(((float)(i+1))/ssidx[l]*4)*4;
		disp[i] = dispF[idx] ;
//		printf("%f %f\n", dispF, disp[i].to_float());
		fprintf(fo, "%.8f %.8f %.8f\n", pyrFFT[i].real().to_float(), pyrFFT[i].imag().to_float(), disp[i].to_float());

	}
	fclose(fo);
	hls::stream<t_input_complex>  sig[C];
	hls::stream<t_disp_scalar>    dispFilter;
	hls::stream<t_nufft_output_complex> sigStreamOut[C];

	fo = fopen("pyr_disp_out.txt", "wb");
	int nRound = 2;
	for(int round=0;round< nRound;round++) {
		for(int i=0;i<ssidx[l];i++) {
			for(int c=0;c<C;c++) {
				sig[c].write(pyrFFT[i]);
			}
			dispFilter.write(disp[i]);
		}
	}
	for(int round=0;round< nRound;round++) {

		nufft_top_t( sig,
				dispFilter,
				sigStreamOut);
	}

	for(int round=0;round< nRound;round++) {

		for(int i=0;i<ssidx[l]*2;i++) {
			t_nufft_output_complex vals[C];
			for(int c=0;c<C;c++) {
				vals[c] = sigStreamOut[c].read();
			}
			fprintf(fo, "%.8f %.8f\n", vals[0].real().to_float(), vals[0].imag().to_float());
		}
	}
	fclose(fo);
}



#endif



#if 0
void nufft_cosfilter1(hls::stream< t_nufft_output_complex > nufftIn[1],
					 hls::stream< t_nufft_output_complex > nufftOut[1],
					 const int nL, const int m) ;
void nufft_filter_tb() {
	hls::stream< t_nufft_output_complex > nufftIn[1];
	hls::stream< t_nufft_output_complex > nufftOut[1];


	FILE * fi=fopen("pyrinBeforeFilterFFT.txt", "rb");
	float vreal,vimag;
	for(int i=0;i<1024;i++) {
		fscanf(fi, "%f %f", &vreal, &vimag);
		nufftIn[0].write( t_nufft_output_complex( vreal, vimag) );
	}
	nufft_cosfilter1(nufftIn, nufftOut, 512, 2);
	for(int i=0;i<512;i++) {
		t_nufft_output_complex val = nufftOut[0].read();
		printf("%.8f %.8f\n", val.real().to_float(), val.imag().to_float());
	}
	fclose(fi);
}
#endif


#include "nufft.h"
//#include "../NUFFT_ForwardFFT/nufft_forwardfft.h"
//#include "pyrrecon.h"

#define C  3
void nufft_top_pyr(hls::stream<t_input_complex>  sig[C],
			   hls::stream<t_disp_scalar>    &dispFilter,
			   hls::stream<t_nufft_output_complex> sigStreamOutH[C],
			   hls::stream<t_nufft_output_complex> sigStreamOutL0[C],
			   hls::stream<t_nufft_output_complex> sigStreamOutLA[C],
			   hls::stream<t_nufft_output_complex> sigStreamOutLP[C]);

void nufft_tb2() {
	float pyr[1520*2];
	float disp[1520];

	int nRound = 4;

	FILE * fi = fopen("NUFFTWARP_PYRIN.dat", "rb");
	fseek(fi, 1520*2 * 59 * sizeof(float), SEEK_SET);
	fread(pyr, 1520*2, sizeof(float), fi);
	fclose(fi);

	fi = fopen("NUFFTWARP_DISIN.dat", "rb");
	fseek(fi, 1520 * 59 * sizeof(float), SEEK_SET);
	fread(disp, 1520, sizeof(float), fi);
	fclose(fi);

	hls::stream<t_input_complex>  sig[C];
	hls::stream<t_disp_scalar>    dispFilter;
	hls::stream<t_nufft_output_complex> sigStreamOutH[C];
	hls::stream<t_nufft_output_complex> sigStreamOutL0[C];
	hls::stream<t_nufft_output_complex> sigStreamOutLA[C];
	hls::stream<t_nufft_output_complex> sigStreamOutLP[C];



//	printf("%.8f %.8f\n", c0[1][361].real().to_float(), c0[1][361].imag().to_float());

	for(int round =0;round<nRound;round++){
		for(int i=0;i<1520;i++) {
			for(int c=0;c<C;c++) {
				t_input_complex v;
				if (i<1504)
					v = t_input_complex( pyr[i*2+0], pyr[i*2+1]);
				else {
					v = t_input_complex( pyr[i*2+0]/64, pyr[i*2+1]/64);
				}
				sig[c].write(v);
			}
			dispFilter.write(disp[i]);
		}
	}
	for(int round =0;round<nRound;round++) {
		nufft_top_pyr(sig, dispFilter, sigStreamOutH, sigStreamOutL0, sigStreamOutLA, sigStreamOutLP);
//		nufft_stream(true, sigStreamOutH, sigStreamOutL0, sigStreamOutLA,
//						   sigStreamFFTOutH, sigStreamFFTOutL0, sigStreamFFTOutLA);
	}


	FILE *fo = fopen("pyr_disp_all_out.txt", "wb");
	for(int round =0;round<nRound;round++) {
		for(int c=0;c<C;c++) {
			for(int i=0;i<1024;i++) {
				t_nufft_output_complex v = sigStreamOutH[c].read();
				if (c == 0) {
	//				std::cout << v << std::endl;
					fprintf(fo, "%.8f %.8f\n", v.real().to_float(), v.imag().to_float());
				}
			}
			fprintf(fo, "Start Low Pass Level 1");
			for(int i=0;i<1024;i++) {
				t_nufft_output_complex v = 	sigStreamOutL0[c].read();
				if (c == 0) {
					//std::cout << v << std::endl;
					fprintf(fo, "%.8f %.8f\n", v.real().to_float(), v.imag().to_float());
				}
			}
			fprintf(fo, "Start Low Pass Level 2-6");
			for(int i=0;i<960;i++) {
				t_nufft_output_complex v = 	sigStreamOutLA[c].read();
				if (c == 0) {
					//std::cout << v << std::endl;
					fprintf(fo, "%.8f %.8f\n", v.real().to_float(), v.imag().to_float());
				}
			}
			fprintf(fo, "Start Low Pass LowerBound");
			for(int i=0;i<16;i++) {
				t_nufft_output_complex v = 	sigStreamOutLP[c].read();
				if (c == 0) {
					//std::cout << v << std::endl;
					fprintf(fo, "%.8f %.8f\n", v.real().to_float(), v.imag().to_float());
				}
			}	}
	}
	fclose(fo);
	printf("Done\n");
}

#if 0

void nufft_full_tb2() {
	float pyr[1520*2];
	float disp[1520];

	int nRound = 4;

	FILE * fi = fopen("NUFFTWARP_PYRIN.dat", "rb");
	fseek(fi, 1520*2 * 59 * sizeof(float), SEEK_SET);
	fread(pyr, 1520*2, sizeof(float), fi);
	fclose(fi);

	fi = fopen("NUFFTWARP_DISIN.dat", "rb");
	fseek(fi, 1520 * 59 * sizeof(float), SEEK_SET);
	fread(disp, 1520, sizeof(float), fi);
	fclose(fi);

	hls::stream<t_input_complex>  sig[C];
	hls::stream<t_disp_scalar>    dispFilter;
	hls::stream<t_nufft_output_complex> sigStreamOutH[C];
	hls::stream<t_nufft_output_complex> sigStreamOutL0[C];
	hls::stream<t_nufft_output_complex> sigStreamOutLA[C];
	hls::stream<t_nufft_output_complex> sigStreamOutLP[C];


	hls::stream<t_nufft_output_complex> sigStreamHOutFFT[C];
	hls::stream<t_nufft_output_complex> sigStreamL0OutFFT[C];
	hls::stream<t_nufft_output_complex> sigStreamLAOutFFT[C];

	hls::stream<t_nufft_output_complex> sigStreamHOutCOSFilter[C];
	hls::stream<t_nufft_output_complex> sigStreamL0OutCOSFilter[C];
	hls::stream<t_nufft_output_complex> sigStreamLAOutCOSFilter[C];

	hls::stream< cmpxDataOut>           sigPyrCOSFilterIn[C];
	hls::stream< t_output_scalar>       sigPyrOut[C];

//	printf("%.8f %.8f\n", c0[1][361].real().to_float(), c0[1][361].imag().to_float());

	for(int round =0;round<nRound;round++){
		for(int i=0;i<1520;i++) {
			for(int c=0;c<C;c++) {
				t_input_complex v;
				if (i<1504)
					v = t_input_complex( pyr[i*2+0], pyr[i*2+1]);
				else {
					v = t_input_complex( pyr[i*2+0]/64, pyr[i*2+1]/64);
				}
				sig[c].write(v);
			}
			dispFilter.write(disp[i]);
		}
	}
	for(int round =0;round<nRound;round++) {
		nufft_top_pyr(sig, dispFilter, sigStreamOutH, sigStreamOutL0, sigStreamOutLA, sigStreamOutLP);

		nufft_stream(true, sigStreamOutH, sigStreamOutL0, sigStreamOutLA,
				sigStreamHOutFFT, sigStreamL0OutFFT, sigStreamLAOutFFT);

		nufft_cosfilter<3>(sigStreamHOutFFT, sigStreamHOutCOSFilter, 512, 2);
		nufft_cosfilter<3>(sigStreamL0OutFFT, sigStreamL0OutCOSFilter, 512, 2);

		for(int k=0;k<4;k++) {
			int sk = 256 >> k;
			nufft_cosfilter<3>(sigStreamLAOutFFT, sigStreamLAOutCOSFilter, sk, 2);
		}
		for(int i=0;i<1520;i++) {
			for(int c=0;c<3;c++) {
				t_nufft_output_complex val;
				cmpxDataOut  valOut;
				if (i<512) {
					val = sigStreamHOutCOSFilter[c].read();
					valOut = cmpxDataOut( val.real(), val.imag());
				}
				if (i>=512 && i < 1024) {
					val = sigStreamL0OutCOSFilter[c].read();
					valOut = cmpxDataOut( val.real(), val.imag());
				}
				if (i>=1024 && i < 1504) {
					val = sigStreamLAOutCOSFilter[c].read();
					valOut = cmpxDataOut( val.real(), val.imag());
				}
				if (i>=1504) {
					val = sigStreamOutLP[c].read();
					valOut = cmpxDataOut( val.real()<<6, val.imag()<<6);
//					valOut = cmpxDataOut( val.real(), val.imag());
				}
				sigPyrCOSFilterIn[c].write(valOut);
			}
		}
		pyr_recon_combine(sigPyrCOSFilterIn, sigPyrOut);

	}

	FILE *fo = fopen("pyr_disp_all_out_recon.txt", "wb");
	for(int round =0;round<nRound;round++) {
		for(int i=0;i<512;i++) {
			for(int c=0;c<3;c++) {
				t_output_scalar val = sigPyrOut[c].read();
				if (c==0)
					fprintf(fo, "%.8f\n", val.to_float());
			}
		}
	}

	fclose(fo);
/*
	FILE *fo = fopen("pyr_disp_all_out_fft_cos.txt", "wb");
	for(int round =0;round<nRound;round++) {
		for(int c=0;c<C;c++) {
			for(int i=0;i<512;i++) {
				t_nufft_output_complex v = sigStreamHOutCOSFilter[c].read();
				if (c == 0) {
	//				std::cout << v << std::endl;
					fprintf(fo, "%.8f %.8f\n", v.real().to_float(), v.imag().to_float());
				}
			}
			for(int i=0;i<512;i++) {
				t_nufft_output_complex v = 	sigStreamL0OutCOSFilter[c].read();
				if (c == 0) {
					//std::cout << v << std::endl;
					fprintf(fo, "%.8f %.8f\n", v.real().to_float(), v.imag().to_float());
				}
			}

			for(int i=0;i<480;i++) {
				t_nufft_output_complex v = 	sigStreamLAOutCOSFilter[c].read();
				if (c == 0) {
					//std::cout << v << std::endl;
					fprintf(fo, "%.8f %.8f\n", v.real().to_float(), v.imag().to_float());
				}
			}
			for(int i=0;i<16;i++) {
				t_nufft_output_complex v = 	sigStreamOutLP[c].read();
				if (c == 0) {
					//std::cout << v << std::endl;
					fprintf(fo, "%.8f %.8f\n", v.real().to_float(), v.imag().to_float());
				}

			}
		}
	}
	fclose(fo);
*/
}
#endif 0

int main() {
//	nufft_tb3();
//	nufft_tb();
	//nufft_filter_tb();
//	pyrrecon_tb();
	nufft_tb2();
//	nufft_full_tb2();
	return 0;
}
