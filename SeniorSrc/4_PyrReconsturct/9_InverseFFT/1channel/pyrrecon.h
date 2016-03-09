#ifndef __PYRRECON__
#define __PYRRECON__

#include "../common_type.h"
#include <hls_fft.h>
//#include "../PyramidConStream/pyr.h"


typedef ap_fixed<16,1> data_in_t;
typedef ap_fixed<24, 24 - 9 + 1 > data_out_t;

typedef std::complex<data_in_t> cmpxDataIn;
typedef std::complex<data_out_t> cmpxDataOut;

typedef ap_fixed<24, 24 - 14 + 1 > data_out_t2;
typedef std::complex<data_out_t2> cmpxDataOut2;


typedef ap_fixed<24,1>            ifft_in_t;
typedef ap_fixed<32, 32 - 22 + 1 > ifft_out_t;

typedef std::complex<ifft_in_t>  cmpxIFFTIn;
typedef std::complex<ifft_out_t> cmpxIFFTOut;


#define  reconC    3

struct config2_recon : hls::ip_fft::params_t {
    static const unsigned ordering_opt = hls::ip_fft::natural_order;
    static const unsigned config_width = 8;

    static const unsigned input_width = 22;
    static const unsigned output_width = 32;

    static const unsigned scaling_opt = hls::ip_fft::unscaled;

    static const unsigned max_nfft = 9;
    static const bool has_nfft = false;
    static const unsigned  channels = 1;

    static const unsigned arch_opt = hls::ip_fft::pipelined_streaming_io;
    //static const unsigned arch_opt = hls::ip_fft::radix_2_burst_io;
};

typedef hls::ip_fft::config_t<config2_recon> config2_recon_t;
typedef hls::ip_fft::status_t<config2_recon> status2_recon_t;


static const data_out_t reconsFilters[1520] = {   0.0000000e+00,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   3.5269635e-02,   7.0256100e-02,   1.0495323e-01,   1.3935523e-01,   1.7345670e-01,   2.0725258e-01,   2.4073816e-01,   2.7390905e-01,   3.0676117e-01,   3.3929075e-01,   3.7149432e-01,   4.0336867e-01,   4.3491085e-01,   4.6611819e-01,   4.9698826e-01,   5.2751883e-01,   5.5770794e-01,   5.8755383e-01,   6.1705493e-01,   6.4620988e-01,   6.7501752e-01,   7.0347686e-01,   7.3158708e-01,   7.5934752e-01,   7.8675771e-01,   8.1381731e-01,   8.4052611e-01,   8.6688408e-01,   8.9289128e-01,   9.1854792e-01,   9.4385433e-01,   9.6881096e-01,   9.9341836e-01,   1.0176772e+00,   1.0415882e+00,   1.0651523e+00,   1.0883704e+00,   1.1112435e+00,   1.1337727e+00,   1.1559593e+00,   1.1778046e+00,   1.1993098e+00,   1.2204763e+00,   1.2413058e+00,   1.2617996e+00,   1.2819595e+00,   1.3017870e+00,   1.3212838e+00,   1.3404517e+00,   1.3592926e+00,   1.3778082e+00,   1.3960004e+00,   1.4138712e+00,   1.4314225e+00,   1.4486563e+00,   1.4655746e+00,   1.4821796e+00,   1.4984732e+00,   1.5144576e+00,   1.5301350e+00,   1.5455074e+00,   1.5605772e+00,   1.5753463e+00,   1.5898172e+00,   1.6039920e+00,   1.6178730e+00,   1.6314625e+00,   1.6447626e+00,   1.6577757e+00,   1.6705042e+00,   1.6829502e+00,   1.6951161e+00,   1.7070043e+00,   1.7186171e+00,   1.7299567e+00,   1.7410255e+00,   1.7518259e+00,   1.7623602e+00,   1.7726307e+00,   1.7826398e+00,   1.7923897e+00,   1.8018829e+00,   1.8111217e+00,   1.8201084e+00,   1.8288453e+00,   1.8373347e+00,   1.8455790e+00,   1.8535805e+00,   1.8613414e+00,   1.8688642e+00,   1.8761510e+00,   1.8832042e+00,   1.8900261e+00,   1.8966188e+00,   1.9029848e+00,   1.9091261e+00,   1.9150452e+00,   1.9207442e+00,   1.9262253e+00,   1.9314907e+00,   1.9365427e+00,   1.9413835e+00,   1.9460152e+00,   1.9504401e+00,   1.9546601e+00,   1.9586777e+00,   1.9624947e+00,   1.9661135e+00,   1.9695360e+00,   1.9727644e+00,   1.9758008e+00,   1.9786473e+00,   1.9813059e+00,   1.9837786e+00,   1.9860676e+00,   1.9881748e+00,   1.9901022e+00,   1.9918518e+00,   1.9934257e+00,   1.9948257e+00,   1.9960538e+00,   1.9971121e+00,   1.9980023e+00,   1.9987265e+00,   1.9992864e+00,   1.9996841e+00,   1.9999213e+00,   0.0000000e+00,   1.9999213e+00,   1.9996841e+00,   1.9992864e+00,   1.9987265e+00,   1.9980023e+00,   1.9971121e+00,   1.9960538e+00,   1.9948257e+00,   1.9934257e+00,   1.9918518e+00,   1.9901022e+00,   1.9881748e+00,   1.9860676e+00,   1.9837786e+00,   1.9813059e+00,   1.9786473e+00,   1.9758008e+00,   1.9727644e+00,   1.9695360e+00,   1.9661135e+00,   1.9624947e+00,   1.9586777e+00,   1.9546601e+00,   1.9504401e+00,   1.9460152e+00,   1.9413835e+00,   1.9365427e+00,   1.9314907e+00,   1.9262253e+00,   1.9207442e+00,   1.9150452e+00,   1.9091261e+00,   1.9029848e+00,   1.8966188e+00,   1.8900261e+00,   1.8832042e+00,   1.8761510e+00,   1.8688642e+00,   1.8613414e+00,   1.8535805e+00,   1.8455790e+00,   1.8373347e+00,   1.8288453e+00,   1.8201084e+00,   1.8111217e+00,   1.8018829e+00,   1.7923897e+00,   1.7826398e+00,   1.7726307e+00,   1.7623602e+00,   1.7518259e+00,   1.7410255e+00,   1.7299567e+00,   1.7186171e+00,   1.7070043e+00,   1.6951161e+00,   1.6829502e+00,   1.6705042e+00,   1.6577757e+00,   1.6447626e+00,   1.6314625e+00,   1.6178730e+00,   1.6039920e+00,   1.5898172e+00,   1.5753463e+00,   1.5605772e+00,   1.5455074e+00,   1.5301350e+00,   1.5144576e+00,   1.4984732e+00,   1.4821796e+00,   1.4655746e+00,   1.4486563e+00,   1.4314225e+00,   1.4138712e+00,   1.3960004e+00,   1.3778082e+00,   1.3592926e+00,   1.3404517e+00,   1.3212838e+00,   1.3017870e+00,   1.2819595e+00,   1.2617996e+00,   1.2413058e+00,   1.2204763e+00,   1.1993098e+00,   1.1778046e+00,   1.1559593e+00,   1.1337727e+00,   1.1112435e+00,   1.0883704e+00,   1.0651523e+00,   1.0415882e+00,   1.0176772e+00,   9.9341836e-01,   9.6881096e-01,   9.4385433e-01,   9.1854792e-01,   8.9289128e-01,   8.6688408e-01,   8.4052611e-01,   8.1381731e-01,   7.8675771e-01,   7.5934752e-01,   7.3158708e-01,   7.0347686e-01,   6.7501752e-01,   6.4620988e-01,   6.1705493e-01,   5.8755383e-01,   5.5770794e-01,   5.2751883e-01,   4.9698826e-01,   4.6611819e-01,   4.3491085e-01,   4.0336867e-01,   3.7149432e-01,   3.3929075e-01,   3.0676117e-01,   2.7390905e-01,   2.4073816e-01,   2.0725258e-01,   1.7345670e-01,   1.3935523e-01,   1.0495323e-01,   7.0256100e-02,   3.5269635e-02,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   0.0000000e+00,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   7.0256100e-02,   1.3935523e-01,   2.0725258e-01,   2.7390905e-01,   3.3929075e-01,   4.0336867e-01,   4.6611819e-01,   5.2751883e-01,   5.8755383e-01,   6.4620988e-01,   7.0347686e-01,   7.5934752e-01,   8.1381731e-01,   8.6688408e-01,   9.1854792e-01,   9.6881096e-01,   1.0176772e+00,   1.0651523e+00,   1.1112435e+00,   1.1559593e+00,   1.1993098e+00,   1.2413058e+00,   1.2819595e+00,   1.3212838e+00,   1.3592926e+00,   1.3960004e+00,   1.4314225e+00,   1.4655746e+00,   1.4984732e+00,   1.5301350e+00,   1.5605772e+00,   1.5898172e+00,   1.6178730e+00,   1.6447626e+00,   1.6705042e+00,   1.6951161e+00,   1.7186171e+00,   1.7410255e+00,   1.7623602e+00,   1.7826398e+00,   1.8018829e+00,   1.8201084e+00,   1.8373347e+00,   1.8535805e+00,   1.8688642e+00,   1.8832042e+00,   1.8966188e+00,   1.9091261e+00,   1.9207442e+00,   1.9314907e+00,   1.9413835e+00,   1.9504401e+00,   1.9586777e+00,   1.9661135e+00,   1.9727644e+00,   1.9786473e+00,   1.9837786e+00,   1.9881748e+00,   1.9918518e+00,   1.9948257e+00,   1.9971121e+00,   1.9987265e+00,   1.9996841e+00,   2.0000000e+00,   1.9996890e+00,   1.9987656e+00,   1.9972443e+00,   1.9951391e+00,   1.9924640e+00,   1.9892326e+00,   1.9854584e+00,   1.9811547e+00,   1.9763344e+00,   1.9710103e+00,   1.9651951e+00,   1.9589012e+00,   1.9521405e+00,   1.9449251e+00,   1.9372668e+00,   1.9291770e+00,   1.9206670e+00,   1.9117480e+00,   1.9024309e+00,   1.8927263e+00,   1.8826448e+00,   1.8721967e+00,   1.8613920e+00,   1.8502409e+00,   1.8387529e+00,   1.8269377e+00,   1.8148046e+00,   1.8023629e+00,   1.7896215e+00,   1.7765894e+00,   1.7632751e+00,   1.7496872e+00,   1.7358341e+00,   1.7217239e+00,   1.7073646e+00,   1.6927642e+00,   1.6779303e+00,   1.6628704e+00,   1.6475920e+00,   1.6321023e+00,   1.6164085e+00,   1.6005174e+00,   1.5844360e+00,   1.5681709e+00,   1.5517286e+00,   1.5351156e+00,   1.5183381e+00,   1.5014024e+00,   1.4843144e+00,   1.4670800e+00,   1.4497050e+00,   1.4321951e+00,   1.4145559e+00,   1.3967927e+00,   1.3789108e+00,   1.3609155e+00,   1.3428119e+00,   1.3246049e+00,   1.3062994e+00,   1.2879002e+00,   1.2694120e+00,   1.2508393e+00,   1.2321866e+00,   1.2134584e+00,   1.1946587e+00,   1.1757920e+00,   1.1568622e+00,   1.1378735e+00,   1.1188296e+00,   1.0997344e+00,   1.0805918e+00,   1.0614053e+00,   1.0421786e+00,   1.0229151e+00,   1.0036184e+00,   9.8429168e-01,   9.6493831e-01,   9.4556149e-01,   9.2616434e-01,   9.0674993e-01,   8.8732126e-01,   8.6788125e-01,   8.4843276e-01,   8.2897859e-01,   8.0952146e-01,   7.9006404e-01,   7.7060894e-01,   7.5115870e-01,   7.3171581e-01,   7.1228268e-01,   6.9286170e-01,   6.7345517e-01,   6.5406535e-01,   6.3469445e-01,   6.1534462e-01,   5.9601796e-01,   5.7671651e-01,   5.5744229e-01,   5.3819723e-01,   5.1898324e-01,   4.9980217e-01,   4.8065584e-01,   4.6154602e-01,   4.4247441e-01,   4.2344270e-01,   4.0445253e-01,   3.8550548e-01,   3.6660311e-01,   3.4774694e-01,   3.2893842e-01,   3.1017901e-01,   2.9147009e-01,   2.7281303e-01,   2.5420915e-01,   2.3565974e-01,   2.1716605e-01,   1.9872930e-01,   1.8035069e-01,   1.6203135e-01,   1.4377242e-01,   1.2557498e-01,   1.0744010e-01,   8.9368795e-02,   7.1362074e-02,   5.3420907e-02,   3.5546233e-02,   1.7738969e-02,   0.0000000e+00,   1.7738969e-02,   3.5546233e-02,   5.3420907e-02,   7.1362074e-02,   8.9368795e-02,   1.0744010e-01,   1.2557498e-01,   1.4377242e-01,   1.6203135e-01,   1.8035069e-01,   1.9872930e-01,   2.1716605e-01,   2.3565974e-01,   2.5420915e-01,   2.7281303e-01,   2.9147009e-01,   3.1017901e-01,   3.2893842e-01,   3.4774694e-01,   3.6660311e-01,   3.8550548e-01,   4.0445253e-01,   4.2344270e-01,   4.4247441e-01,   4.6154602e-01,   4.8065584e-01,   4.9980217e-01,   5.1898324e-01,   5.3819723e-01,   5.5744229e-01,   5.7671651e-01,   5.9601796e-01,   6.1534462e-01,   6.3469445e-01,   6.5406535e-01,   6.7345517e-01,   6.9286170e-01,   7.1228268e-01,   7.3171581e-01,   7.5115870e-01,   7.7060894e-01,   7.9006404e-01,   8.0952146e-01,   8.2897859e-01,   8.4843276e-01,   8.6788125e-01,   8.8732126e-01,   9.0674993e-01,   9.2616434e-01,   9.4556149e-01,   9.6493831e-01,   9.8429168e-01,   1.0036184e+00,   1.0229151e+00,   1.0421786e+00,   1.0614053e+00,   1.0805918e+00,   1.0997344e+00,   1.1188296e+00,   1.1378735e+00,   1.1568622e+00,   1.1757920e+00,   1.1946587e+00,   1.2134584e+00,   1.2321866e+00,   1.2508393e+00,   1.2694120e+00,   1.2879002e+00,   1.3062994e+00,   1.3246049e+00,   1.3428119e+00,   1.3609155e+00,   1.3789108e+00,   1.3967927e+00,   1.4145559e+00,   1.4321951e+00,   1.4497050e+00,   1.4670800e+00,   1.4843144e+00,   1.5014024e+00,   1.5183381e+00,   1.5351156e+00,   1.5517286e+00,   1.5681709e+00,   1.5844360e+00,   1.6005174e+00,   1.6164085e+00,   1.6321023e+00,   1.6475920e+00,   1.6628704e+00,   1.6779303e+00,   1.6927642e+00,   1.7073646e+00,   1.7217239e+00,   1.7358341e+00,   1.7496872e+00,   1.7632751e+00,   1.7765894e+00,   1.7896215e+00,   1.8023629e+00,   1.8148046e+00,   1.8269377e+00,   1.8387529e+00,   1.8502409e+00,   1.8613920e+00,   1.8721967e+00,   1.8826448e+00,   1.8927263e+00,   1.9024309e+00,   1.9117480e+00,   1.9206670e+00,   1.9291770e+00,   1.9372668e+00,   1.9449251e+00,   1.9521405e+00,   1.9589012e+00,   1.9651951e+00,   1.9710103e+00,   1.9763344e+00,   1.9811547e+00,   1.9854584e+00,   1.9892326e+00,   1.9924640e+00,   1.9951391e+00,   1.9972443e+00,   1.9987656e+00,   1.9996890e+00,   2.0000000e+00,   1.9996841e+00,   1.9987265e+00,   1.9971121e+00,   1.9948257e+00,   1.9918518e+00,   1.9881748e+00,   1.9837786e+00,   1.9786473e+00,   1.9727644e+00,   1.9661135e+00,   1.9586777e+00,   1.9504401e+00,   1.9413835e+00,   1.9314907e+00,   1.9207442e+00,   1.9091261e+00,   1.8966188e+00,   1.8832042e+00,   1.8688642e+00,   1.8535805e+00,   1.8373347e+00,   1.8201084e+00,   1.8018829e+00,   1.7826398e+00,   1.7623602e+00,   1.7410255e+00,   1.7186171e+00,   1.6951161e+00,   1.6705042e+00,   1.6447626e+00,   1.6178730e+00,   1.5898172e+00,   1.5605772e+00,   1.5301350e+00,   1.4984732e+00,   1.4655746e+00,   1.4314225e+00,   1.3960004e+00,   1.3592926e+00,   1.3212838e+00,   1.2819595e+00,   1.2413058e+00,   1.1993098e+00,   1.1559593e+00,   1.1112435e+00,   1.0651523e+00,   1.0176772e+00,   9.6881096e-01,   9.1854792e-01,   8.6688408e-01,   8.1381731e-01,   7.5934752e-01,   7.0347686e-01,   6.4620988e-01,   5.8755383e-01,   5.2751883e-01,   4.6611819e-01,   4.0336867e-01,   3.3929075e-01,   2.7390905e-01,   2.0725258e-01,   1.3935523e-01,   7.0256100e-02,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   0.0000000e+00,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.3935523e-01,   2.7390905e-01,   4.0336867e-01,   5.2751883e-01,   6.4620988e-01,   7.5934752e-01,   8.6688408e-01,   9.6881096e-01,   1.0651523e+00,   1.1559593e+00,   1.2413058e+00,   1.3212838e+00,   1.3960004e+00,   1.4655746e+00,   1.5301350e+00,   1.5898172e+00,   1.6447626e+00,   1.6951161e+00,   1.7410255e+00,   1.7826398e+00,   1.8201084e+00,   1.8535805e+00,   1.8832042e+00,   1.9091261e+00,   1.9314907e+00,   1.9504401e+00,   1.9661135e+00,   1.9786473e+00,   1.9881748e+00,   1.9948257e+00,   1.9987265e+00,   2.0000000e+00,   1.9987656e+00,   1.9951391e+00,   1.9892326e+00,   1.9811547e+00,   1.9710103e+00,   1.9589012e+00,   1.9449251e+00,   1.9291770e+00,   1.9117480e+00,   1.8927263e+00,   1.8721967e+00,   1.8502409e+00,   1.8269377e+00,   1.8023629e+00,   1.7765894e+00,   1.7496872e+00,   1.7217239e+00,   1.6927642e+00,   1.6628704e+00,   1.6321023e+00,   1.6005174e+00,   1.5681709e+00,   1.5351156e+00,   1.5014024e+00,   1.4670800e+00,   1.4321951e+00,   1.3967927e+00,   1.3609155e+00,   1.3246049e+00,   1.2879002e+00,   1.2508393e+00,   1.2134584e+00,   1.1757920e+00,   1.1378735e+00,   1.0997344e+00,   1.0614053e+00,   1.0229151e+00,   9.8429168e-01,   9.4556149e-01,   9.0674993e-01,   8.6788125e-01,   8.2897859e-01,   7.9006404e-01,   7.5115870e-01,   7.1228268e-01,   6.7345517e-01,   6.3469445e-01,   5.9601796e-01,   5.5744229e-01,   5.1898324e-01,   4.8065584e-01,   4.4247441e-01,   4.0445253e-01,   3.6660311e-01,   3.2893842e-01,   2.9147009e-01,   2.5420915e-01,   2.1716605e-01,   1.8035069e-01,   1.4377242e-01,   1.0744010e-01,   7.1362074e-02,   3.5546233e-02,   0.0000000e+00,   3.5546233e-02,   7.1362074e-02,   1.0744010e-01,   1.4377242e-01,   1.8035069e-01,   2.1716605e-01,   2.5420915e-01,   2.9147009e-01,   3.2893842e-01,   3.6660311e-01,   4.0445253e-01,   4.4247441e-01,   4.8065584e-01,   5.1898324e-01,   5.5744229e-01,   5.9601796e-01,   6.3469445e-01,   6.7345517e-01,   7.1228268e-01,   7.5115870e-01,   7.9006404e-01,   8.2897859e-01,   8.6788125e-01,   9.0674993e-01,   9.4556149e-01,   9.8429168e-01,   1.0229151e+00,   1.0614053e+00,   1.0997344e+00,   1.1378735e+00,   1.1757920e+00,   1.2134584e+00,   1.2508393e+00,   1.2879002e+00,   1.3246049e+00,   1.3609155e+00,   1.3967927e+00,   1.4321951e+00,   1.4670800e+00,   1.5014024e+00,   1.5351156e+00,   1.5681709e+00,   1.6005174e+00,   1.6321023e+00,   1.6628704e+00,   1.6927642e+00,   1.7217239e+00,   1.7496872e+00,   1.7765894e+00,   1.8023629e+00,   1.8269377e+00,   1.8502409e+00,   1.8721967e+00,   1.8927263e+00,   1.9117480e+00,   1.9291770e+00,   1.9449251e+00,   1.9589012e+00,   1.9710103e+00,   1.9811547e+00,   1.9892326e+00,   1.9951391e+00,   1.9987656e+00,   2.0000000e+00,   1.9987265e+00,   1.9948257e+00,   1.9881748e+00,   1.9786473e+00,   1.9661135e+00,   1.9504401e+00,   1.9314907e+00,   1.9091261e+00,   1.8832042e+00,   1.8535805e+00,   1.8201084e+00,   1.7826398e+00,   1.7410255e+00,   1.6951161e+00,   1.6447626e+00,   1.5898172e+00,   1.5301350e+00,   1.4655746e+00,   1.3960004e+00,   1.3212838e+00,   1.2413058e+00,   1.1559593e+00,   1.0651523e+00,   9.6881096e-01,   8.6688408e-01,   7.5934752e-01,   6.4620988e-01,   5.2751883e-01,   4.0336867e-01,   2.7390905e-01,   1.3935523e-01,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   0.0000000e+00,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   2.7390905e-01,   5.2751883e-01,   7.5934752e-01,   9.6881096e-01,   1.1559593e+00,   1.3212838e+00,   1.4655746e+00,   1.5898172e+00,   1.6951161e+00,   1.7826398e+00,   1.8535805e+00,   1.9091261e+00,   1.9504401e+00,   1.9786473e+00,   1.9948257e+00,   2.0000000e+00,   1.9951391e+00,   1.9811547e+00,   1.9589012e+00,   1.9291770e+00,   1.8927263e+00,   1.8502409e+00,   1.8023629e+00,   1.7496872e+00,   1.6927642e+00,   1.6321023e+00,   1.5681709e+00,   1.5014024e+00,   1.4321951e+00,   1.3609155e+00,   1.2879002e+00,   1.2134584e+00,   1.1378735e+00,   1.0614053e+00,   9.8429168e-01,   9.0674993e-01,   8.2897859e-01,   7.5115870e-01,   6.7345517e-01,   5.9601796e-01,   5.1898324e-01,   4.4247441e-01,   3.6660311e-01,   2.9147009e-01,   2.1716605e-01,   1.4377242e-01,   7.1362074e-02,   0.0000000e+00,   7.1362074e-02,   1.4377242e-01,   2.1716605e-01,   2.9147009e-01,   3.6660311e-01,   4.4247441e-01,   5.1898324e-01,   5.9601796e-01,   6.7345517e-01,   7.5115870e-01,   8.2897859e-01,   9.0674993e-01,   9.8429168e-01,   1.0614053e+00,   1.1378735e+00,   1.2134584e+00,   1.2879002e+00,   1.3609155e+00,   1.4321951e+00,   1.5014024e+00,   1.5681709e+00,   1.6321023e+00,   1.6927642e+00,   1.7496872e+00,   1.8023629e+00,   1.8502409e+00,   1.8927263e+00,   1.9291770e+00,   1.9589012e+00,   1.9811547e+00,   1.9951391e+00,   2.0000000e+00,   1.9948257e+00,   1.9786473e+00,   1.9504401e+00,   1.9091261e+00,   1.8535805e+00,   1.7826398e+00,   1.6951161e+00,   1.5898172e+00,   1.4655746e+00,   1.3212838e+00,   1.1559593e+00,   9.6881096e-01,   7.5934752e-01,   5.2751883e-01,   2.7390905e-01,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   0.0000000e+00,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   5.2751883e-01,   9.6881096e-01,   1.3212838e+00,   1.5898172e+00,   1.7826398e+00,   1.9091261e+00,   1.9786473e+00,   2.0000000e+00,   1.9811547e+00,   1.9291770e+00,   1.8502409e+00,   1.7496872e+00,   1.6321023e+00,   1.5014024e+00,   1.3609155e+00,   1.2134584e+00,   1.0614053e+00,   9.0674993e-01,   7.5115870e-01,   5.9601796e-01,   4.4247441e-01,   2.9147009e-01,   1.4377242e-01,   0.0000000e+00,   1.4377242e-01,   2.9147009e-01,   4.4247441e-01,   5.9601796e-01,   7.5115870e-01,   9.0674993e-01,   1.0614053e+00,   1.2134584e+00,   1.3609155e+00,   1.5014024e+00,   1.6321023e+00,   1.7496872e+00,   1.8502409e+00,   1.9291770e+00,   1.9811547e+00,   2.0000000e+00,   1.9786473e+00,   1.9091261e+00,   1.7826398e+00,   1.5898172e+00,   1.3212838e+00,   9.6881096e-01,   5.2751883e-01,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   0.0000000e+00,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   9.6881096e-01,   1.5898172e+00,   1.9091261e+00,   2.0000000e+00,   1.9291770e+00,   1.7496872e+00,   1.5014024e+00,   1.2134584e+00,   9.0674993e-01,   5.9601796e-01,   2.9147009e-01,   0.0000000e+00,   2.9147009e-01,   5.9601796e-01,   9.0674993e-01,   1.2134584e+00,   1.5014024e+00,   1.7496872e+00,   1.9291770e+00,   2.0000000e+00,   1.9091261e+00,   1.5898172e+00,   9.6881096e-01,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.2246468e-16,   1.0000000e+00,   1.0000000e+00,   1.0000000e+00,   1.0000000e+00,   1.0000000e+00,   8.7484360e-01,   6.0672918e-01,   2.9800898e-01,   0.0000000e+00,   2.9800898e-01,   6.0672918e-01,   8.7484360e-01,   1.0000000e+00,   1.0000000e+00,   1.0000000e+00,   1.0000000e+00 };

void pyr_recon(hls::stream< cmpxDataOut>  &pyrFFT,
		       hls::stream< t_recon_complex>  &imDFTOut);

void pyr_recon_ifft(hls::stream< t_recon_complex>  &DFTOut,  hls::stream<  PIXEL_RAW>  &RGBOut);

//void pyr_recon_combine(hls::stream< cmpxDataOut>  pyrFFT[reconC],
//				  hls::stream< PIXEL_RAW>  &sigOut);

#endif
