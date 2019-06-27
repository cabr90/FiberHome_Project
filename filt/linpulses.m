function y = linpulses(a,Nxs,pls)
%LINPULSES genera delle sequenze di impulsi modulati linearmente 
% 
%  y=LINPULSES(a,Nxs,pls) genera N segnali campionati con Nxs punti per
%  tempo di simbolo, ciascuno costituito da un treno di impulsi di forma data e
%  modulati linearmente secondo la sequenza a:
%
%    y[k]=sum_i(a[i]*p(k*T/Nxs-i*T),  
% 
%  dove k il tempo discreto e p(t) è l'impulso scelto, opportunamente normalizzato per avere
%  energia unitaria: int(|p(t)|^2)=int(|P(f)|^2)=1
%  Si assume che la trasformata dell'impulso sia nulla al di fuori della
%  banda di frequenze [-Nxs/2T,Nxs/2T].
%
%   *** ARGOMENTI DELLA FUNZIONE ***
%  a:       sequenza di simboli modulanti (righe=modi, colonne=simboli)
%  Nxs:     numero di campioni per tempo di simbolo
%  pls:     caratteristiche dell'impulso (pls.shape,pls.ord,pls.bw)

T          = 1;
dt         = T/Nxs;
[m,Na]     = size(a);
df         = 1/(T*Na);
N          = Na*Nxs;
N1         = floor(N/2);
N2         = ceil(N/2)-1;
f          = df*[0:N2,-N1:-1];   %vettore delle frequenze dopo la FFT

H=trans_func(pls,f);
%if (N2<N1)
%    H(N1+1)=0.;
%end

%normalizzazione per avere energia unitaria
norm = sqrt(N/dt/sum(abs(H).^2));
H    = H*norm;

%sequenza di impulsi di tipo delta (con Nxs campioni per simbolo) modulati
%secondo i simboli dati
y            = complex(zeros(m,N));
y(:,1:Nxs:N) = a;

%shaping spettrale secondo l'impulso desiderato
Y       = fft(y,[],2);
Y       = Y.*H;
y       = ifft(Y,[],2);

% H_match = conj(H);
% u       =  ifft(fft(y,[],2).*H_match);
% y_rx    =  u(:,1:Nxs:end);

end