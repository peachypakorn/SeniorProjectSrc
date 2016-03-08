
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include "forwardfft512.h"
#include <stdio.h>

#include <fstream>
#include <string>
#include <sstream>

using namespace std;
// Temporary Declaration
typedef ap_fixed<17, 4> t_nufft_output_scalar;
typedef std::complex< t_nufft_output_scalar > t_nufft_output_complex;


#define BUF_SIZE 64
int main()
{
	static hls::stream< t_nufft_output_complex > xn_input;
		static hls::stream< t_nufft_output_complex> xk_output;
		FILE * fi = fopen("inputMem.txt", " wb");
int ndat = 512 + 256+128+64;//960
//		int ndat = 256 ;
		for(int j=0;j<4;j++)
		{
			for(int i=0;i<ndat;i++) {
				t_nufft_output_complex val ((float)( rand()^0xFFFF) /65535, (float)( rand()^0xFFFF) /65535);
				if(j==0){
					fprintf(fi,"%.8f %.8f \n",val.real().to_float(),val.imag().to_float());
						}
				xn_input.write( val);
				//xn_input[1].write( val);
				//xn_input[2].write( val);
			}
			//s2_fft_512_stream3( xn_input[0], xn_input[1], xn_input[2],
				//	xk_output[0], xk_output[1], xk_output[2]);+
			s2_fft_512_stream(xn_input,xk_output);
			FILE * fo = fopen("pyr_out.txt", " wb");
			for(int i=0;i<ndat/2;i++) {
				t_nufft_output_complex valOut1 = xk_output.read();
				//t_nufft_output_complex valOut2 =xk_output[1].read();
				//t_nufft_output_complex valOut3 =xk_output[2].read();
				fprintf(fo, " %.8f %.8f  ", valOut1.real().to_double(), valOut1.imag().to_double());
				//fprintf(fo, "%.8f %.8f  ", valOut2.real().to_double(), valOut2.imag().to_double());
				//fprintf(fo, "%.8f %.8f \n", valOut3.real().to_double(), valOut3.imag().to_double());

				std::cout << valOut1/* << valOut2 << valOut3*/<< std::endl;
			}
			fclose(fo);
		}
		fclose(fi);
					return 0;
};
