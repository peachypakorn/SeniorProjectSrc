#include "../common_type.h"

template < int NUMOUT>
void  S2_S3DeMultiplex( 	hls::stream<t_nufft_output_complex > S2Out_H[3],
							hls::stream<t_nufft_output_complex > S2Out_L0[3],
							hls::stream<t_nufft_output_complex > S2Out_LA[3],
							hls::stream<t_nufft_output_complex > S2Out_LP[3],

							hls::stream<t_nufft_output_complex > S3In_H[NUMOUT][3],
							hls::stream<t_nufft_output_complex > S3In_L0[NUMOUT][3],
							hls::stream<t_nufft_output_complex > S3In_LA[NUMOUT][3],
							hls::stream<t_nufft_output_complex > S3In_LP[NUMOUT][3]) {
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
			S3In_H[inputC][0].write( S2Out_H[0].read());
			S3In_H[inputC][1].write( S2Out_H[1].read());
			S3In_H[inputC][2].write( S2Out_H[2].read());
		}
	}

	for(int inputC=0;inputC<NUMOUT;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<512;i++) {
#pragma HLS PIPELINE
			S3In_L0[inputC][0].write( S2Out_L0[0].read());
			S3In_L0[inputC][1].write( S2Out_L0[1].read());
			S3In_L0[inputC][2].write( S2Out_L0[2].read());
		}
	}

	for(int inputC=0;inputC<NUMOUT;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<256+128+64+32;i++) {
#pragma HLS PIPELINE
			S3In_LA[inputC][0].write( S2Out_LA[0].read());
			S3In_LA[inputC][1].write( S2Out_LA[1].read());
			S3In_LA[inputC][2].write( S2Out_LA[2].read());
		}
	}


	for(int inputC=0;inputC<NUMOUT;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<16;i++) {
#pragma HLS PIPELINE
			S3In_LP[inputC][0].write( S2Out_LP[0].read());
			S3In_LP[inputC][1].write( S2Out_LP[1].read());
			S3In_LP[inputC][2].write( S2Out_LP[2].read());
		}
	}

}


void  S2_S3DeMultiplex_2( hls::stream<t_nufft_output_complex > S2Out_H[3],
		hls::stream<t_nufft_output_complex > S2Out_L0[3],
		hls::stream<t_nufft_output_complex > S2Out_LA[3],
		hls::stream<t_nufft_output_complex > S2Out_LP[3],

		hls::stream<t_nufft_output_complex > S3In_H[2][3],
		hls::stream<t_nufft_output_complex > S3In_L0[2][3],
		hls::stream<t_nufft_output_complex > S3In_LA[2][3],
		hls::stream<t_nufft_output_complex > S3In_LP[2][3]) {

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
