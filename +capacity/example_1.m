clear all;

M        = 64;                      % MQAM
SNRdB    = (-10:10).';                % Es/N0 [dB]
SNR      = 10.^(SNRdB/10);

rate     = [4];                  % overall rate (H(P))

stars    = qammod((0:M-1)',M,'UnitAveragePower',true,'PlotConstellation',false);


% compute the uniform distribution for the MQAM constellation %
Pu       = ones(M,1)./M;

% compute the Maxwel-Bolzan distribution of the MQAM for the specified
% entropy
PMB      = zeros(M,1);
lmd      = fzero(@(x) entropy(maxbolz(stars,stars,x))-rate,[0,6]);
PMB(:,1) = maxbolz(stars,stars,lmd);

% calculate the symbols and the noise avg energies
Es_u     = capacity_functions.symbol_energy(stars,Pu);
sg_u     = sqrt(Es_u) .* 10.^(-SNRdB/20);

Es_mb    = capacity_functions.symbol_energy(stars,PMB);
sg_mb    = sqrt(Es_mb) .* 10.^(-SNRdB/20);

% estimates the gmi for the MQAM symbols through the Gauss-Hermite
% quadrature method
gmi_ps_u     = arrayfun(@(x) capacity_functions.qam_gmi(stars,M,x,Pu),sg_u);
gmi_ps_pmb   = arrayfun(@(x) capacity_functions.qam_gmi(stars,M,x,PMB),sg_mb);

% estimates the mi for the MQAM symbols through the Gauss-Hermite
% quadrature method
mi_ps_u      = arrayfun(@(x) capacity_functions.qam_mi(stars,M,x,Pu),sg_u);
mi_ps_pmb    = arrayfun(@(x) capacity_functions.qam_mi(stars,M,x,PMB),sg_mb);

% estimates the gmi and the mi for the MQAM symbols through the Monte-Carlo
% integration
D          = 10000;                    
for i = 1:length(SNRdB)
    z  = wgn(D,1,sg_mb(i)^2,'linear','complex');
    zu = wgn(D,1,sg_u(i)^2,'linear','complex');
    gmi_ps_u_mc(i,1)   = capacity_functions.qam_montecarlo_gmi(zu,stars,M,sg_u(i),Pu);
    mi_ps_u_mc(i,1)    = capacity_functions.qam_montecarlo_mi(zu,stars,M,sg_u(i),Pu);
    gmi_ps_pmb_mc(i,1) = capacity_functions.qam_montecarlo_gmi(z,stars,M,sg_mb(i),PMB);
    mi_ps_pmb_mc(i,1)  = capacity_functions.qam_montecarlo_mi(z,stars,M,sg_mb(i),PMB);
end

figure(1);
subplot(2,2,1); grid on; hold on; box on;
title('GMI - PS-MQAM');
plot(SNRdB,gmi_ps_u,SNRdB,gmi_ps_pmb,SNRdB,gmi_ps_pmb_mc)
xlabel('E_s/N_0 [dB]');
ylabel('GMI [bit/QAM symbol]')
legend('64-QAM uniform','PS-64-QAM with GH','PS-64-QAM with MC')

subplot(2,2,2); grid on; hold on; box on;
title('MI - PS-MQAM');
plot(SNRdB,mi_ps_u,SNRdB,mi_ps_pmb,SNRdB,mi_ps_pmb_mc)
xlabel('E_s/N_0 [dB]');
ylabel('MI [bit/QAM symbol]')
legend('64-QAM uniform','PS-64-QAM with GH','PS-64-QAM with MC')

subplot(2,2,3); grid on; hold on; box on;
title('MI - GMI with GH - PS-MQAM');
plot(SNRdB,(mi_ps_u-gmi_ps_u),SNRdB,(mi_ps_pmb-gmi_ps_pmb))
xlabel('E_s/N_0 [dB]');
ylabel('MI - GMI [bit/QAM symbol]')
legend('MI-GMI uniform-64-QAM with GH','MI-GMI PS-64-QAM with GH')

subplot(2,2,4); grid on; hold on; box on;
title('MI - GMI with MC - PS-MQAM');
plot(SNRdB,mi_ps_u_mc-gmi_ps_u_mc,SNRdB,mi_ps_pmb_mc-gmi_ps_pmb_mc)
xlabel('E_s/N_0 [dB]');
ylabel('MI - GMI [bit/QAM symbol]')
legend('MI-GMI uniform-64-QAM with MC','MI-GMI PS-64-QAM with MC')
