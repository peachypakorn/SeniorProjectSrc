#ifndef __FORWARDFFT1024__
#define __FORWARDFFT1024__

#include "ap_fixed.h"
#include "hls_fft.h"

#include "../common_type.h"

#include <complex>


// Temporary Declaration
//typedef ap_fixed<17, 4> t_nufft_output_scalar;
//typedef std::complex< t_nufft_output_scalar > t_nufft_output_complex;



void s2_fft_1024_stream3(
		hls::stream< t_nufft_output_complex > &in,
		hls::stream< t_nufft_output_complex > &out) ;


#endif
