function y = match_samp(u,Nxs,pls)

T           = 1;
dt          = T/Nxs;
[n,Nu]      = size(u);
df          = 1/(T*Nu/Nxs);       % Spaziatura tra componenti frequenziali (normalizzata alla symbol rate)
N           = Nu;
N1          = floor(N/2);
N2          = ceil(N/2)-1;
f           = df*[0:N2,-N1:-1];   %vettore delle frequenze dopo la FFT



Hf_match    = conj(trans_func(pls,f));

% norm        = sqrt(N*dt/sum(abs(Hf_match).^2));
% Hf_match    = Hf_match*norm;

u           =  ifft(fft(u,[],2).*Hf_match,[],2);
y           =  u(:,1:Nxs:end);
end