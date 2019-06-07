clear all;

addpath filt;
import capacity.*;

lqam                    = [2];


sig.Rs                  = 40e9;                  % symbol rate
sig.N_simboli           = 2^18;                  % Lunghezza della sequenza di simboli trasmessi
sig.Nxs                 = 2;                     % Discrete-time representation: samples per symbol
sig.filt_tx.shape       = 'RRC';                 % Transmission filter shape (transfer function)
sig.filt_tx.rate        = 1;                     % Symbol rate for which the TX filter satisfies Nyquist condition (normalized to actual symbol rate Rs)
sig.filt_tx.roff        = 0.2;                   % Roll-off factor of TX filter
sig.QAM_order           = 2^lqam(1);             % Inizializza il numero di livelli QAM

wss.filt_wss.shape      = 'WSS';
wss.filt_wss.BWotf      = 8.5e9;
wss.filt_wss.B          = 50e9;
wss.filt_wss.rate       = sig.Rs;
wss.N_wss               = 40;

SNRdB                   = (-5:2:30).';           % Es/N0 [dB]
SNR                     = 10.^(SNRdB/10);

for k = 1:length(SNRdB)
    stars                  = qammod((0:sig.QAM_order-1)',sig.QAM_order,'UnitAveragePower',true,'PlotConstellation',false);
    
    
    c_tx                   = qammod(randi([0 sig.QAM_order-1],sig.N_simboli,1),sig.QAM_order,'UnitAveragePower',true,'PlotConstellation',false);
    p_qam                  = ones(sig.QAM_order,1)*1/sig.QAM_order;
    
    Es                     = capacity_functions.symbol_energy(stars,p_qam);
    sg                     = sqrt(0.5*sig.Nxs*Es).* 10.^(-SNRdB(k)/20);
    
    %% modulazione + rumore awgn + filtro matchato
    % Modulazione lineare dell'impulso base con rrc roll-off alpha
    u_tx                   = linpulses(c_tx.',sig.Nxs,sig.filt_tx);
    % rumore gaussiona
    z                      = sg*(randn(1,length(u_tx))+1i*randn(1,length(u_tx)));
    u                      = u_tx + z;
    
    %% AIR estimation only match filter
    % Filtro matchato e campionamento
    y_rx                   = match_samp(u,sig.Nxs,sig.filt_tx).';
    % Stima parametri
    h0                     = mean(real(y_rx.*conj(c_tx)))/mean(abs(c_tx).^2);
    sgn                    = sqrt(var(y_rx) - abs(h0)^2*var(c_tx));
    % Stima MI 
    mi_mc(k,1)             = capacity_functions.qam_montecarlo_mi(y_rx-h0*c_tx,h0*stars,sig.QAM_order,sgn,p_qam);
    
    %% AIR estimation WSS + match filter
    % Filtro WSS
    u_wss                  = wss_filter(u,sig.Nxs,sig.Rs,wss.filt_wss,wss.N_wss);
    % Filtro matchato e campionamento a tempo di simbolo
    y_rx_wss               = match_samp(u_wss,sig.Nxs,sig.filt_tx).';
    % Stima parametri
    h0_wss                 = mean(real(y_rx_wss.*conj(c_tx)))/mean(abs(c_tx).^2);
    sgn_wss                = sqrt(var(y_rx_wss) - abs(h0_wss)^2*var(c_tx));
    % Stima MI 
    mi_mc_wss(k,1)         = capacity_functions.qam_montecarlo_mi    (y_rx_wss-h0_wss*c_tx,h0_wss*stars,sig.QAM_order,sgn_wss,p_qam);
    
%     mi_mc_sim_wss(k,1)     = capacity_functions.qam_montecarlo_mi_sim(y_rx_wss,h0_wss*c_tx,h0_wss*stars,sig.QAM_order,sgn_wss,p_qam);
    %% MI awgn gauss hermit estimation with wss
    sg_gh                  = sqrt(Es).* 10.^(-SNRdB(k)/20);
    h0_wss                 = mean(real(y_rx_wss.*conj(c_tx)))/mean(abs(c_tx).^2);
    sgn_wss                = sqrt(var(y_rx_wss) - abs(h0_wss)^2*var(c_tx));
    mi_gh(k,1)             = capacity_functions.qam_mi(h0_wss*stars,sig.QAM_order,sgn_wss,p_qam);
    
end

hold on; grid on; legend;
plot(SNRdB,mi_mc        ,'-sr', 'DisplayName',horzcat('    Nxs = ',num2str(sig.Nxs)));
plot(SNRdB,mi_mc_wss    ,'--ob','DisplayName',horzcat('wss Nxs = ',num2str(sig.Nxs)));
% plot(SNRdB,mi_mc_sim_wss,'--oy','DisplayName',horzcat('wss Nxs = ',num2str(sig.Nxs)));
plot(SNRdB,mi_gh        ,'--*g','DisplayName',horzcat('Gauss-Hermit'));