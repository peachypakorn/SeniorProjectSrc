function writeFloat( filename, pyrIn)
fidd = fopen(filename,'w');
fwrite(fidd, [real(pyrIn(:))],'float32', 'n');

fclose(fidd);