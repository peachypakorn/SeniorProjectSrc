#include "pyrrecon.h"

void S3_ReconMergeOld(hls::stream< t_nufft_output_complex > H[reconC],
					hls::stream< t_nufft_output_complex > L0[reconC],
				hls::stream< t_nufft_output_complex > LA[reconC],
				hls::stream< t_nufft_output_complex > LP[reconC],
		        hls::stream< t_recon_complex>  imDFTOut[reconC]) {
#pragma HLS inline off

#pragma HLS dataflow

#pragma HLS INTERFACE axis port=H
#pragma HLS INTERFACE axis port=L0
#pragma HLS INTERFACE axis port=LA
#pragma HLS INTERFACE axis port=LP
#pragma HLS INTERFACE axis port=imDFTOut

#pragma HLS data_pack variable=H
#pragma HLS data_pack variable=L0
#pragma HLS data_pack variable=LA
#pragma HLS data_pack variable=LP
#pragma HLS data_pack variable=imDFTOut

#pragma HLS STREAM variable=H   depth=1024
#pragma HLS STREAM variable=L0  depth=1024
#pragma HLS STREAM variable=LA  depth=960
#pragma HLS STREAM variable=LP  depth=32

	hls::stream< cmpxDataOut>           sigPyrCOSFilterIn[reconC];


#pragma HLS data_pack variable=sigPyrCOSFilterIn

	for(int i=0;i<1520;i++) {
#pragma HLS pipeline
		for(int c=0;c<reconC;c++) {
			t_nufft_output_complex val;
			if (i<512)
				val = H[c].read();
			if (i>=512 && i < 1024)
				val = L0[c].read();
			if (i>=1024 && i < 1504)
				val = LA[c].read();
			if (i>=1504) {
				val = LP[c].read();
				//t_nufft_output_complex vali = LP[c].read();
				//val.real() = vali.real() << 8;
				//val.imag() = vali.imag() << 8;
			}

			cmpxDataOut  valOut( val.real(), val.imag());
			sigPyrCOSFilterIn[c].write(valOut);
		}
	}
	pyr_recon(sigPyrCOSFilterIn, imDFTOut );
}



#define K 2
void S3_ReconMerge(hls::stream< t_nufft_output_complex > H[reconC],
					hls::stream< t_nufft_output_complex > L0[reconC],
				hls::stream< t_nufft_output_complex > LA[reconC],
				hls::stream< t_nufft_output_complex > LP[reconC],
		        hls::stream< t_recon_complex>  imDFTOut[reconC]) {
#pragma HLS inline off

#pragma HLS dataflow

#pragma HLS INTERFACE axis port=H
#pragma HLS INTERFACE axis port=L0
#pragma HLS INTERFACE axis port=LA
#pragma HLS INTERFACE axis port=LP
#pragma HLS INTERFACE axis port=imDFTOut

#pragma HLS data_pack variable=H
#pragma HLS data_pack variable=L0
#pragma HLS data_pack variable=LA
#pragma HLS data_pack variable=LP
#pragma HLS data_pack variable=imDFTOut

#pragma HLS STREAM variable=H   depth=1024
#pragma HLS STREAM variable=L0  depth=1024
#pragma HLS STREAM variable=LA  depth=960
#pragma HLS STREAM variable=LP  depth=32

	hls::stream< cmpxDataOut>           sigPyrCOSFilterIn[K][reconC];
#pragma HLS data_pack variable=sigPyrCOSFilterIn
#pragma HLS STREAM variable=sigPyrCOSFilterIn depth=1520

	for(int k=0;k<K;k++) {
		for(int i=0;i<1520;i++) {
	#pragma HLS pipeline
			for(int c=0;c<reconC;c++) {
				t_nufft_output_complex val;
				if (i<512)
					val = H[c].read();
				if (i>=512 && i < 1024)
					val = L0[c].read();
				if (i>=1024 && i < 1504)
					val = LA[c].read();
				if (i>=1504) {
					val = LP[c].read();
					//t_nufft_output_complex vali = LP[c].read();
					//val.real() = vali.real() << 8;
					//val.imag() = vali.imag() << 8;
				}

				cmpxDataOut  valOut( val.real(), val.imag());
				sigPyrCOSFilterIn[k][c].write(valOut);
			}
		}
	}
	hls::stream< t_recon_complex>  imDFTOuts[2][reconC];
#pragma HLS data_pack variable=imDFTOuts
#pragma HLS STREAM variable=imDFTOuts depth=1024

	pyr_recon(sigPyrCOSFilterIn[0], imDFTOuts[0] );
	pyr_recon(sigPyrCOSFilterIn[1], imDFTOuts[1] );

	for(int k=0;k<K;k++) {
		for(int i=0;i<512;i++) {
#pragma HLS PIPELINE
			for(int c=0;c<reconC;c++) {
				imDFTOut[c].write( imDFTOuts[k][c].read());
			}
		}
	}

}
