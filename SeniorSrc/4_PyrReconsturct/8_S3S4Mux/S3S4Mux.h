#ifndef __S3S4_MUX
#define __S3S4_MUX

#include "../common_type.h"

void S3S4Mux_2(hls::stream< t_recon_complex>  S3Out[2][3],
			 	 hls::stream< t_recon_complex>  DFTOut[3]);


#endif


