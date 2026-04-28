function U0 = mraf_backward_fft(Uf)
% mraf_backward_fft Inverse of mraf_forward_fft.
N = size(Uf, 1);
U0 = ifftshift(ifft2(ifftshift(Uf))) * N;
end