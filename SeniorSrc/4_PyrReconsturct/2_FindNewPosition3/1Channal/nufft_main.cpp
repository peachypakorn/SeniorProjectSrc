#include "nufft.h"
#include "nufft_func.h"

template<typename T, typename T2>
T cmpconv( const T2 &in) {
#pragma HLS inline
	T2 valIn = in;
	T val(valIn.real(), valIn.imag());
return val;

}

#define C  3

void nufft_top_t(hls::stream<t_input_complex>  &sig,
			   hls::stream<t_disp_scalar>    &dispFilter,
			   hls::stream<t_nufft_output_complex> &sigStreamOut) {
#pragma HLS data_pack variable=sig
#pragma HLS data_pack variable=sigStreamOut



#pragma HLS DATAFLOW
	nufft_top<C, 512>( sig, dispFilter, sigStreamOut, 512,0, 1023);
}

const int Kset = 7;
const int limits[] = { 512, 512,256,128,64,32,16};

const int Klimits[] = { 255, 255, 127, 63, 31, 15, 3};
const int climits[] = { 0, 512, 1024,1280,1408,1472,1504};
const int llimits[] = { 9,9,8,7,6,5,4};
const int lshifts[] = { 1, 2, 3,3,3,5,9};
const int pshifts[] = { 8, 7, 5,4,3,0, -5};



void nufft_top_pyr(hls::stream<t_input_complex>  &sig,
			   hls::stream<t_disp_scalar>    &dispFilter,
			   hls::stream<t_nufft_output_complex> &sigStreamOutH,
			   hls::stream<t_nufft_output_complex> &sigStreamOutL0,
			   hls::stream<t_nufft_output_complex> &sigStreamOutLA,
			   hls::stream<t_nufft_output_complex> &sigStreamOutLP) {
#pragma HLS inline off
/*
#ifndef  NUFFTB
#pragma HLS INTERFACE axis port=sig
#pragma HLS INTERFACE axis port=dispFilter
#pragma HLS INTERFACE axis port=sigStreamOutH
#pragma HLS INTERFACE axis port=sigStreamOutL0
#pragma HLS INTERFACE axis port=sigStreamOutLA
#pragma HLS INTERFACE axis port=sigStreamOutLP
#endif
*/

#pragma HLS data_pack variable=sig

#pragma HLS DATAFLOW
	hls::stream<t_input_complex>  sigH;
	hls::stream<t_input_complex>  sigL0;
	hls::stream<t_input_complex>  sigLA;

	hls::stream<t_disp_scalar > disp0;
	hls::stream<t_disp_scalar > disp1;
	hls::stream<t_disp_scalar > disp2;



#pragma HLS data_pack variable=sigH
#pragma HLS data_pack variable=sigL0
#pragma HLS data_pack variable=sigLA
#pragma HLS data_pack variable=sigStreamOutH
#pragma HLS data_pack variable=sigStreamOutL0
#pragma HLS data_pack variable=sigStreamOutLA
#pragma HLS data_pack variable=sigStreamOutLP

#pragma HLS stream variable=disp0   depth=512
#pragma HLS stream variable=disp1   depth=512
#pragma HLS stream variable=disp2   depth=512

#pragma HLS stream variable=sigH    depth=512
#pragma HLS stream variable=sigL0   depth=512
#pragma HLS stream variable=sigLA   depth=512
//#pragma HLS stream variable=sigStreamOutH    depth=512
//#pragma HLS stream variable=sigStreamOutL0   depth=512
//#pragma HLS stream variable=sigStreamOutLA   depth=490
//#pragma HLS stream variable=sigStreamOutLP   depth=64


	hls::stream<t_input_complex>  sigInMem;
#pragma HLS data_pack variable=sigInMem
#pragma HLS stream variable=sigInMem    depth=1520
	for(int coefIdx = 0;coefIdx < 1520; coefIdx++) {
#pragma HLS pipeline
		sigInMem.write(sig.read());
	}

	int l = 0;
	int i = 0;
	for(int coefIdx = 0;coefIdx < 1520; coefIdx++)
	{
#pragma HLS pipeline

			t_input_complex v = sigInMem.read();
			if (coefIdx >=0 && coefIdx < climits[1]) 			sigH.write(v);
			if (coefIdx >=climits[1] && coefIdx < climits[2])	sigL0.write(v);
			if (coefIdx >=climits[2] && coefIdx < climits[6])   sigLA.write(v);
			if (coefIdx >=climits[6])                           sigStreamOutLP.write(v);


		t_disp_scalar dispVal = dispFilter.read();

		if (coefIdx >=0 && coefIdx < climits[1])   	        disp0.write(dispVal);
		if (coefIdx >=climits[1] && coefIdx < climits[2]) 	disp1.write(dispVal);
		if (coefIdx >=climits[2] && coefIdx < climits[6]) 	disp2.write(dispVal);
	}
part1:
	nufft_top<C, 512>( sigH, disp0,  sigStreamOutH, 512,255);
part2:
	nufft_top<C, 512>( sigL0, disp1, sigStreamOutL0, 512,255);

part3:

//const int limits[] = { 512, 512,256,128,64,32,16};
//const int Klimits[] = { 255, 255, 127, 63, 31, 15, 3};

	for(int k=0;k<4;k++) {
//#pragma HLS DATAFLOW
		const int limit = 256>>k;
		const int klimit = 127 >> k;
		int level = k;
		nufft_top<C, 256>(sigLA, disp2, sigStreamOutLA, limit, level, klimit);
	}

//	for(int k=2;k<6;k++) {
//		nufft_top<C, 256>(sigL[2], disp2, sigStreamOutLA, limits[k],k-1,Klimits[k]);
//	}


	//if (disp0.)
}

