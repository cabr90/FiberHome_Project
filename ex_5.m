clear all;
import capacity.*;
import pshape.*;

% Set parameters
N           = 2^12;             % number of samples
M           = 4;                % modulation order MQAM
SNRdB       = (-5:15).';        % Es/N0 [dB]
SNR         = 10.^(SNRdB/10);
D           = 1;                % Complex dimension (real dimension 2D)

% Probabilistic shaping parameters
make_shape  = false;            % set true if you want to use prob shaping  
block_size  = 2^12;             % size of the shaped block
N_blocks    = N/block_size;     % number of blocks
rate        = 3.98;             % rate of the modulation 
n           = block_size;       % block size
dm          = 'ccdm';           % distribution matcher type (ccdm or ess)

% Construct the shaper and the constellation points
if make_shape
    ps      = pshape(dm,M,n,rate);
end
[stars,p]   = get_stars(M,D);


for i = 1:length(SNRdB)
    
    if make_shape
        [x,bits_tx,p] = get_ps_samples(ps,D,N_blocks);
    else
        x             = qammod(randi([0 M-1],N,D),M,'UnitAveragePower',true);
    end
     
    % calculate the symbols and the noise avg energies
    Es     = capacity_functions.symbol_energy(stars,p);
    sg     = sqrt(1/(2*D)*Es) .* 10.^(-SNRdB/20);
    
    % AWGN channel
    z           = sg(i)*(randn(length(x),D)+1i*randn(length(x),D));
    y           = x + z;
    
    % capacity for a D dimensional complex signal    
    mi_mc (i,1)  = capacity_functions.md_mi_montecarlo(z,stars,M^D,sg(i),p);
    gmi_mc(i,1)  = capacity_functions.md_gmi_montecarlo(z,stars,M^D,sg(i),p);
    
    % capacity for a 2D dimensional real signal
    z2d = []; stars2d = []; p2d = [];
    for j = 1:D
        z2d         = [z2d,real(z(:,j)),imag(z(:,j))];
        stars2d     = [stars2d,real(stars(:,j)),imag(stars(:,j))];
        p2d         = p;
    end
    
    fun1            = @ custom1_get_index_kth_bit;
    mi_mc_2d(i,1)   = capacity_functions.md_mi_montecarlo(z2d,stars2d,M^D,sg(i),p2d);
    gmi_mc_2d(i,1)  = capacity_functions.md_gmi_montecarlo(z2d,stars2d,M^D,sg(i),p2d,'idx_fun',fun1);
end

hold on; grid on; legend;
plot(SNRdB,mi_mc ,'-r');
plot(SNRdB,mi_mc_2d ,'-or');
plot(SNRdB,gmi_mc ,'-b');
plot(SNRdB,gmi_mc_2d ,'-ob');

function [x,b,p]   = get_ps_samples(ps,D,N_blocks)

    c_tx                    = [];
    bits_tx                 = [];
    
    for j = 1:D*N_blocks
        [c_tx_tmp,bits_tx_tmp,px]  = ps.encoding;
        c_tx                       = [c_tx;c_tx_tmp];
        bits_tx                    = [bits_tx;bits_tx_tmp];
    end
    
    x     = reshape(c_tx,[],D);
    b     = reshape(bits_tx,[],D);
    
    p     = P_QAM(x,M);

end

function [stars,p] = get_stars(M,D)

% create a M-QAM constellation
stars    = qammod((0:M-1)'.*ones(M,D),M,'UnitAveragePower',true,'PlotConstellation',false);

% compute the uniform distribution for the M-QAM constellation %
p       = ones(M^D,1)./M^D;

% create the M^D points
if D>1
    Ctmp = [];
    for i = 1:M
        for j = 2:D
            Ctmp = [Ctmp;ones(M,1)*stars(i,1:j-1),stars(:,j:end)];
        end
    end
    stars  = Ctmp;
end

end

function idx       = custom1_get_index_kth_bit(k,M,C,b)
[n,d] = size(C);

m   = sqrt(M^(2/d));
B   = [];
idx = [];

n1 = modnorm(C(:,1),'avpow',1);
n2 = modnorm(pammod((0:m-1),m),'avpow',1);

norm = n1/n2;

for  i = 1:d
    B    = [B,de2bi(pamdemod(C(:,i).'*norm,m,0,'gray'),log2(m),'left-msb')];
end

idx  = find(B(:,k)==b);

end

% function idx = custom2_get_index_kth_bit(k,M,C,b)
% [n,d] = size(C);
% 
% m   = M^(2/d);
% B   = [];
% idx = [];
% 
% Ctmp = [];
% for i = 1:2:d
%     Ctmp = [Ctmp,C(:,i)+1i*C(:,i+1)];
% end
% C = Ctmp;
% 
% for  i = 1:d/2
%     B    = [B,reshape(qamdemod(C(:,i).',m,'OutputType','bit','UnitAveragePower',true).',[],log2(m))];
% end
% 
% idx  = find(B(:,k)==b);
% 
% end









