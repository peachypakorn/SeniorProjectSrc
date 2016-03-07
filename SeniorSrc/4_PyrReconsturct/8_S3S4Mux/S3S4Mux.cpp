#include "../common_type.h"
// I won't do a template now.
template < int NUMIN>
void S3S4Mux(hls::stream< t_recon_complex>  S3Out[NUMIN][3],
			 hls::stream< t_recon_complex>  DFTOut[3]) {

#pragma HLS DATAFLOW
#pragma HLS INLINE

#pragma HLS data_pack variable=S3Out
#pragma HLS data_pack variable=DFTOut

	for(int inputC=0;inputC<NUMIN;inputC++) {
//#pragma HLS DATAFLOW
		for(int i=0;i<512;i++) {
#pragma HLS PIPELINE
			//t_nufft_output_complex   T0;
			//t_nufft_output_complex   T1;
			//t_nufft_output_complex   T2;
			DFTOut[0].write( S3Out[inputC][0].read());
			DFTOut[1].write( S3Out[inputC][1].read());
			DFTOut[2].write( S3Out[inputC][2].read());
		}
	}


}
void S3S4Mux_2 (hls::stream< t_recon_complex>  S3Out[2][3],
			 	 hls::stream< t_recon_complex>  DFTOut[3]) {

#pragma HLS DATAFLOW
#pragma HLS INLINE

#pragma HLS INTERFACE axis port=S3Out
#pragma HLS INTERFACE axis port=DFTOut

#pragma HLS data_pack variable=S3Out
#pragma HLS data_pack variable=DFTOut

#pragma HLS stream variable=S3Out  depth=512
#pragma HLS stream variable=DFTOut depth=512

	S3S4Mux<2>(S3Out, DFTOut);

}
