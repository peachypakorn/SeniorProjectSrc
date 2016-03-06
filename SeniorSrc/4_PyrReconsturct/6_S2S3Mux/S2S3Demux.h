#ifndef __S2S3__DEMUX
#define __S2S3__DEMUX
#include "../common_type.h"

void  S2_S3DeMultiplex_2( hls::stream<t_nufft_output_complex > S2Out_H[3],
		hls::stream<t_nufft_output_complex > S2Out_L0[3],
		hls::stream<t_nufft_output_complex > S2Out_LP[3],
		hls::stream<t_nufft_output_complex > S2Out_LA[3],

		hls::stream<t_nufft_output_complex > S3In_H[2][3],
		hls::stream<t_nufft_output_complex > S3In_L0[2][3],
		hls::stream<t_nufft_output_complex > S3In_LP[2][3],
		hls::stream<t_nufft_output_complex > S3In_LA[2][3]);


void  S1_S2Multiplex_4( 	hls::stream<t_nufft_output_complex > S1Out_H[4][3],
						hls::stream<t_nufft_output_complex > S1Out_L0[4][3],
						hls::stream<t_nufft_output_complex > S1Out_LP[4][3],
						hls::stream<t_nufft_output_complex > S1Out_LA[4][3],

						hls::stream<t_nufft_output_complex > S2In_H[3],
						hls::stream<t_nufft_output_complex > S2In_L0[3],
						hls::stream<t_nufft_output_complex > S2In_LP[3],
						hls::stream<t_nufft_output_complex > S2In_LA[3]);

#endif
