#include "../common_type.h"

template < int NUMIN>
void  S1_S2Multiplex( 	hls::stream<t_nufft_output_complex > S1Out_H[NUMIN][3],
						hls::stream<t_nufft_output_complex > S1Out_L0[NUMIN][3],
						hls::stream<t_nufft_output_complex > S1Out_LA[NUMIN][3],
						hls::stream<t_nufft_output_complex > S1Out_LP[NUMIN][3],

						hls::stream<t_nufft_output_complex > S2In_H[3],
						hls::stream<t_nufft_output_complex > S2In_L0[3],
						hls::stream<t_nufft_output_complex > S2In_LA[3],
						hls::stream<t_nufft_output_complex > S2In_LP[3]) {
#pragma HLS DATAFLOW
#pragma HLS INLINE

#pragma HLS data_pack variable=S1Out_H
#pragma HLS data_pack variable=S1Out_L0
#pragma HLS data_pack variable=S1Out_L0
#pragma HLS data_pack variable=S1Out_LP
#pragma HLS data_pack variable=S1Out_LA


#pragma HLS data_pack variable=S2In_H
#pragma HLS data_pack variable=S2In_L0
#pragma HLS data_pack variable=S2In_L0
#pragma HLS data_pack variable=S2In_LP
#pragma HLS data_pack variable=S2In_LA


	for(int inputC=0;inputC<NUMIN;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<1024;i++) {
#pragma HLS PIPELINE
			//t_nufft_output_complex   T0;
			//t_nufft_output_complex   T1;
			//t_nufft_output_complex   T2;
			S2In_H[0].write( S1Out_H[inputC][0].read());
			S2In_H[1].write( S1Out_H[inputC][1].read());
			S2In_H[2].write( S1Out_H[inputC][2].read());
		}
	}


	for(int inputC=0;inputC<NUMIN;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<1024;i++) {
#pragma HLS PIPELINE
			S2In_L0[0].write( S1Out_L0[inputC][0].read());
			S2In_L0[1].write( S1Out_L0[inputC][1].read());
			S2In_L0[2].write( S1Out_L0[inputC][2].read());
		}
	}

	for(int inputC=0;inputC<NUMIN;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<512+256+128+64;i++) {
#pragma HLS PIPELINE
			S2In_LA[0].write( S1Out_LA[inputC][0].read());
			S2In_LA[1].write( S1Out_LA[inputC][1].read());
			S2In_LA[2].write( S1Out_LA[inputC][2].read());
		}
	}


	for(int inputC=0;inputC<NUMIN;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<16;i++) {
#pragma HLS PIPELINE
			S2In_LP[0].write( S1Out_LP[inputC][0].read());
			S2In_LP[1].write( S1Out_LP[inputC][1].read());
			S2In_LP[2].write( S1Out_LP[inputC][2].read());
		}
	}

}


void  S1_S2Multiplex_4( 	hls::stream<t_nufft_output_complex > S1Out_H[4][3],
						hls::stream<t_nufft_output_complex > S1Out_L0[4][3],
						hls::stream<t_nufft_output_complex > S1Out_LA[4][3],
						hls::stream<t_nufft_output_complex > S1Out_LP[4][3],

						hls::stream<t_nufft_output_complex > S2In_H[3],
						hls::stream<t_nufft_output_complex > S2In_L0[3],
						hls::stream<t_nufft_output_complex > S2In_LA[3],
						hls::stream<t_nufft_output_complex > S2In_LP[3]) {

#pragma HLS DATAFLOW

#ifndef  NUFFTB

#pragma HLS INTERFACE axis port=S1Out_H
#pragma HLS INTERFACE axis port=S1Out_L0
#pragma HLS INTERFACE axis port=S1Out_L0
#pragma HLS INTERFACE axis port=S1Out_LP
#pragma HLS INTERFACE axis port=S1Out_LA


#pragma HLS INTERFACE axis port=S2In_H
#pragma HLS INTERFACE axis port=S2In_L0
#pragma HLS INTERFACE axis port=S2In_L0
#pragma HLS INTERFACE axis port=S2In_LP
#pragma HLS INTERFACE axis port=S2In_LA

#endif
#pragma HLS data_pack variable=S1Out_H
#pragma HLS data_pack variable=S1Out_L0
#pragma HLS data_pack variable=S1Out_L0
#pragma HLS data_pack variable=S1Out_LP
#pragma HLS data_pack variable=S1Out_LA


#pragma HLS data_pack variable=S2In_H
#pragma HLS data_pack variable=S2In_L0
#pragma HLS data_pack variable=S2In_LP
#pragma HLS data_pack variable=S2In_LA


#pragma HLS STREAM variable=S2In_H depth=4096 dim=2
#pragma HLS stream variable=S2In_L0 depth=4096 dim=2
#pragma HLS stream variable=S2In_LP depth=4096 dim=2
#pragma HLS stream variable=S2In_LA depth=256 dim=2


 S1_S2Multiplex<4>(S1Out_H, S1Out_L0, S1Out_LA, S1Out_LP,
		 S2In_H, S2In_L0, S2In_LA, S2In_LP);

}
