function writef(filename, variable),
fid = fopen(filename, 'wb');
fprintf(fid, '%3.8f %3.8f\n', real(variable), imag(variable));
fclose(fid);
end