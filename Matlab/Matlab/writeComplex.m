function writeComplex( filename, pyrIn)
fidd = fopen(filename,'w');
fwrite(fidd, [real(pyrIn(:))'; imag(pyrIn(:))'],'float32', 'n');

fclose(fidd);