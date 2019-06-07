function u = wss_filter(u,Nxs,Rs,pls,N_wss)

T      = 1;
dt     = T/Nxs;
Ntot   = length(u);
N1     = floor(Ntot/2);
N2     = ceil(Ntot/2)-1;
df     = Nxs/Ntot;              % Spaziatura tra componenti frequenziali (normalizzata alla symbol rate)
f      = df*[0:N2,-N1:-1]*Rs;   % Vettore delle frequenze per la FFT (normalizzato alla symbol rate)

H_wss  = trans_func(pls,f);     % Funzione di trasferimento di un filtro WSS

H_wss  = H_wss.^(N_wss);

u      = fft(u);                % Trasforma nel dominio della frequenza
u      = H_wss.*u;
u      = ifft(u);

end