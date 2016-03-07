#ifndef __S3__RECON
#define __S3__RECON
#include "pyrrecon.h"

void S3_ReconMerge(hls::stream< t_nufft_output_complex > H[reconC],
					hls::stream< t_nufft_output_complex > L0[reconC],
				hls::stream< t_nufft_output_complex > LA[reconC],
				hls::stream< t_nufft_output_complex > LP[reconC],
		        hls::stream< t_recon_complex>  imDFTOut[reconC]);

#endif
