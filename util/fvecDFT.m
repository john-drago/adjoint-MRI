function f = fvecDFT(fs, N, shift)
% This function will generate a corresponding vector of frequencies for the
% DFT that is performed by MATLAB functions. 
arguments
    fs
    N
    shift = 1
end
% also consider 0:1/N:1-1/N
f_dft = linspace(0,1,N+1);
f_dft = f_dft(1:end-1);
f_dft(f_dft >= 0.5) = f_dft(f_dft >= 0.5) - 1;
f = f_dft * fs;

if shift
    f = fftshift(f);
end

end