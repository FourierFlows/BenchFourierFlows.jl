function qht = RHS_eq(qh)

global invKsq K L Ksq nu nnu

psih = -invKsq.*qh;
uh = -1i*L.*psih;
vh = +1i*K.*psih;

u = real(ifft2(uh));
v = real(ifft2(vh));
q = real(ifft2(qh));

qht = -1i*K.*(fft2(u.*q)) - 1i*L.*(fft2(v.*q)) - nu*Ksq.^nnu.*qh;
end