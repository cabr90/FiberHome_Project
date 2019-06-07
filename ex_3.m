clear all;

addpath filt;
import capacity.*;
import pshape.*;
import pshape.dm.*;


lqam                    = [4];
Nxs                     = [4];


sig.Rs                  = 34e9;                  % symbol rate
sig.N_simboli           = 2^12;                  % Lunghezza della sequenza di simboli trasmessi
sig.Nxs                 = Nxs(1);                % Discrete-time representation: samples per symbol
sig.filt_tx.shape       = 'RRC';                 % Transmission filter shape (transfer function)
sig.filt_tx.rate        = 1;                     % Symbol rate for which the TX filter satisfies Nyquist condition (normalized to actual symbol rate Rs)
sig.filt_tx.roff        = 0.2;                   % Roll-off factor of TX filter
sig.QAM_order           = 2^lqam(1);             % Inizializza il numero di livelli QAM


wss.filt_wss.shape      = 'WSS';
wss.filt_wss.BWotf      = 8.5e9;
wss.filt_wss.B          = 35e9;
wss.filt_wss.rate       = sig.Rs;
wss.N_wss               = 1;

SNRdB                   = (-5:2:30).';           % Es/N0 [dB]
SNR                     = 10.^(SNRdB/10);



block_size           = 2^12;
N_blocks             = sig.N_simboli/block_size;

rate = 2.8;            % rate of the modulation 
N    = block_size;     % block size
M    = sig.QAM_order;  % modulation order
dm   = 'ccdm';         % distribution matcher type (only ccdm now)

ps                      = pshape(dm,M,N,rate);

for k = 1:length(SNRdB)
    
    
    c_tx                    = [];
    bits_tx                 = [];

    
    for i = 1:N_blocks
        % Probabilistic Shaping encoding
        [c_tx_tmp,bits_tx_tmp,px]  = ps.encoding;
        c_tx                    = [c_tx;c_tx_tmp];
        bits_tx                 = [bits_tx;bits_tx_tmp];
    end
    
    stars                  = qammod((0:sig.QAM_order-1)',sig.QAM_order,'UnitAveragePower',true,'PlotConstellation',false);
    p_qam                  = P_QAM(c_tx,sig.QAM_order);
    
%     c_tx                   = qammod(randi([0 sig.QAM_order-1],sig.N_simboli,1),sig.QAM_order,'UnitAveragePower',true,'PlotConstellation',false);
%     p_qam                  = ones(sig.QAM_order,1)*1/sig.QAM_order;
    
    Es                     = capacity_functions.symbol_energy(stars,p_qam);
    sg                     = sqrt(0.5*sig.Nxs*Es).* 10.^(-SNRdB(k)/20);
    
    %% modulazione + rumore awgn + filtro matchato
    % Modulazione lineare dell'impulso base con rrc roll-off alpha
    u_tx                   = linpulses(c_tx.',sig.Nxs,sig.filt_tx);
    % rumore gaussiona
    z                      = sg*(randn(1,length(u_tx))+1i*randn(1,length(u_tx)));
    u                      = u_tx + z;
    % Filtro WSS
%     u                      = wss_filter(u,sig.Nxs,sig.Rs,wss.filt_wss,wss.N_wss);
    % Filtro matchato e campionamento a tempo di simbolo
    y_rx                   = match_samp( u ,sig.Nxs,sig.filt_tx).';
    
    %% AIR estimation with montecarlo integration
    h0                     = real(mean(y_rx.*conj(c_tx))/mean(abs(c_tx).^2));
    sgn                    = sqrt(mean(abs(y_rx).^2) - abs(h0)^2*mean(abs(c_tx).^2));
    mi_mc(k,1)             = capacity_functions.qam_montecarlo_mi(y_rx-h0*c_tx,h0*stars,sig.QAM_order,sgn,p_qam);
    
end

hold on; grid on; legend;
plot(SNRdB,mi_mc,'DisplayName',horzcat('Nxs = ',num2str(sig.Nxs)));
% plot(SNRdB,mi_mc_b2bu,'DisplayName',horzcat('mc uniform Nxs = ',num2str(sig.Nxs))');