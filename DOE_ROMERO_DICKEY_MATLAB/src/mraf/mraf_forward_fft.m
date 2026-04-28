function Uf = mraf_forward_fft(U0)
% mraf_forward_fft Orthogonal ideal FFT propagation used by MRAF.
N = size(U0, 1);
Uf = fftshift(fft2(fftshift(U0))) / N;
end