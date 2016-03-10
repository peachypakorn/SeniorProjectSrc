
// Includes
#include <stdio.h>
#include "cplx_data.h"

// Public functions
void cplx_data_get_string(char* c, cplx_data_t data)
{
//    int whole, thousandths;
//    int whole1, thousandths1;
//    whole =  data.data_re;
//	    thousandths = ( data.data_re - whole) * 1000;
//	    if (thousandths<0) thousandths *=-1;
//	whole1 =  data.data_im;
//	    thousandths1 = ( data.data_im - whole1) * 1000;
//	    if (thousandths1<0) thousandths1 *=-1;
//	    //xil_printf("%d.%03d\n", whole, thousandths);

	//sprintf(c, "%d.%03d + j*%d.%03d", whole, thousandths, whole1, thousandths1);
	//sprintf(c, "%d + j*%d", data.data_re, data.data_im);
	sprintf(c, "%d", data.data_im);
}
