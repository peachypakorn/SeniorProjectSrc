m = 2;
nL = 512;

mySig = fft(f1.nufftBefore);
mySig1 = mySig./ (cos(pi*(fftshift(0:m*nL-1)-m*nL/2)/(m*nL))');
myF = [ mySig1(1:256); mySig1(769:1024)];

mySig2 = [ mySig(1:256) ./ cos(pi * ((0:255)/(m*nL)))';  mySig(769:1024) ./ cos(pi * ((-256:-1)/(m*nL)))' ];
%myF = [ mySig(1:256); mySig(769:1024)];


V = [ cos(pi * ((0:255)/(m*nL)))' ; cos(pi * ((-256:-1)/(m*nL)))'];
Vinv = 1./V;

fftSig = fftshift(fft(f1.nufftBefore));
%fftSig = fftshift((1:m*nL)');
fftSig = fftSig ./ cos(pi*((0:m*nL-1)-m*nL/2)/(m*nL))';
%Fall = fftSig;
F = ifftshift(fftSig(m*nL/2-floor(nL/2)+(1:nL)));
%F = ifft(ifftshift(fftSig(m*nL/2-floor(nL/2)+(1:nL))));