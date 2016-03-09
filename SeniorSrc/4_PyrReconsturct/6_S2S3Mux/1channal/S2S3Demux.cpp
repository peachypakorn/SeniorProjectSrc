#include "../common_type.h"

template < int NUMOUT>
void  S2_S3DeMultiplex( 	hls::stream<t_nufft_output_complex > &S2Out_H,
							hls::stream<t_nufft_output_complex > &S2Out_L0,
							hls::stream<t_nufft_output_complex > &S2Out_LA,
							hls::stream<t_nufft_output_complex > &S2Out_LP,

							hls::stream<t_nufft_output_complex > S3In_H[NUMOUT],
							hls::stream<t_nufft_output_complex > S3In_L0[NUMOUT],
							hls::stream<t_nufft_output_complex > S3In_LA[NUMOUT],
							hls::stream<t_nufft_output_complex > S3In_LP[NUMOUT]) {
#pragma HLS DATAFLOW
#pragma HLS INLINE

#pragma HLS data_pack variable=S2Out_H
#pragma HLS data_pack variable=S2Out_L0
#pragma HLS data_pack variable=S2Out_LP
#pragma HLS data_pack variable=S2Out_LA


#pragma HLS data_pack variable=S3In_H
#pragma HLS data_pack variable=S3In_L0
#pragma HLS data_pack variable=S3In_LP
#pragma HLS data_pack variable=S3In_LA

	for(int inputC=0;inputC<NUMOUT;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<512;i++) {
#pragma HLS PIPELINE
			//t_nufft_output_complex   T0;
			//t_nufft_output_complex   T1;
			//t_nufft_output_complex   T2;
			S3In_H[inputC].write( S2Out_H.read());

		}
	}

	for(int inputC=0;inputC<NUMOUT;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<512;i++) {
#pragma HLS PIPELINE
			S3In_L0[inputC].write( S2Out_L0.read());
			}
	}

	for(int inputC=0;inputC<NUMOUT;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<256+128+64+32;i++) {
#pragma HLS PIPELINE
			S3In_LA[inputC].write( S2Out_LA.read());
		}
	}


	for(int inputC=0;inputC<NUMOUT;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<16;i++) {
#pragma HLS PIPELINE
			S3In_LP[inputC].write( S2Out_LP.read());

		}
	}

}


void  S2_S3DeMultiplex_2( hls::stream<t_nufft_output_complex > &S2Out_H,
		hls::stream<t_nufft_output_complex > &S2Out_L0,
		hls::stream<t_nufft_output_complex > &S2Out_LA,
		hls::stream<t_nufft_output_complex > &S2Out_LP,

		hls::stream<t_nufft_output_complex > S3In_H[2],
		hls::stream<t_nufft_output_complex > S3In_L0[2],
		hls::stream<t_nufft_output_complex > S3In_LA[2],
		hls::stream<t_nufft_output_complex > S3In_LP[2]) {

#pragma HLS DATAFLOW

#pragma HLS INTERFACE axis port=S2Out_H
#pragma HLS INTERFACE axis port=S2Out_L0
#pragma HLS INTERFACE axis port=S2Out_LP
#pragma HLS INTERFACE axis port=S2Out_LA

#pragma HLS stream variable=S2Out_H  depth=512
#pragma HLS stream variable=S2Out_L0 depth=512
#pragma HLS stream variable=S2Out_LP depth=512
#pragma HLS stream variable=S2Out_LA depth=16


#pragma HLS stream variable=S3In_H  depth=512
#pragma HLS stream variable=S3In_L0 depth=512
#pragma HLS stream variable=S3In_LP depth=512
#pragma HLS stream variable=S3In_LA depth=16

#pragma HLS INTERFACE axis port=S3In_H
#pragma HLS INTERFACE axis port=S3In_L0
#pragma HLS INTERFACE axis port=S3In_LP
#pragma HLS INTERFACE axis port=S3In_LA

#pragma HLS data_pack variable=S2Out_H
#pragma HLS data_pack variable=S2Out_L0
#pragma HLS data_pack variable=S2Out_LP
#pragma HLS data_pack variable=S2Out_LA


#pragma HLS data_pack variable=S3In_H
#pragma HLS data_pack variable=S3In_L0
#pragma HLS data_pack variable=S3In_LP
#pragma HLS data_pack variable=S3In_LA


	S2_S3DeMultiplex<2>(S2Out_H, S2Out_L0, S2Out_LA, S2Out_LP,
 		               S3In_H, S3In_L0, S3In_LA, S3In_LP);

}
