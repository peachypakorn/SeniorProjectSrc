#ifndef __S2S3__DEMUX
#define __S2S3__DEMUX
#include "../common_type.h"

void  S2_S3DeMultiplex_2( hls::stream<t_nufft_output_complex > &S2Out_H ,
		hls::stream<t_nufft_output_complex > &S2Out_L0 ,
		hls::stream<t_nufft_output_complex > &S2Out_LP ,
		hls::stream<t_nufft_output_complex > &S2Out_LA ,

		hls::stream<t_nufft_output_complex > S3In_H[2] ,
		hls::stream<t_nufft_output_complex > S3In_L0[2] ,
		hls::stream<t_nufft_output_complex > S3In_LP[2] ,
		hls::stream<t_nufft_output_complex > S3In_LA[2] );


void  S1_S2Multiplex_4( 	hls::stream<t_nufft_output_complex > S1Out_H[4] ,
						hls::stream<t_nufft_output_complex > S1Out_L0[4] ,
						hls::stream<t_nufft_output_complex > S1Out_LP[4] ,
						hls::stream<t_nufft_output_complex > S1Out_LA[4] ,

						hls::stream<t_nufft_output_complex > &S2In_H ,
						hls::stream<t_nufft_output_complex > &S2In_L0 ,
						hls::stream<t_nufft_output_complex > &S2In_LP ,
						hls::stream<t_nufft_output_complex > &S2In_LA );

#endif
