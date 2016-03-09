#include "../S3_ReconMerge/pyrrecon.h"



void S4_finalrecon(hls::stream< t_recon_complex>  DFTOut[reconC],
		hls::stream< PIXEL_RAW>  &sigOut) {
#pragma HLS INTERFACE axis port=DFTOut

	// Enable this if instantiated this alone
#pragma HLS INTERFACE axis port=sigOut

#pragma HLS data_pack variable=DFTOut
#pragma HLS data_pack variable=sigOut

#pragma HLS DATAFLOW

	pyr_recon_ifft(DFTOut, sigOut);

}
