import pshape.*;
%% initialization
rate = 3.58;           % rate of the modulation 
N    = 128;            % block size
M    = 16;             % modulation order
dm   = 'ess';          % distribution matcher type (ccdm or ess)

%% construction of the probabilistic shaping class
ps          = pshape(dm,M,N,rate);
%% encoding
[s,b,p_qam] = ps.encoding('bits_real',randi([0 1], [ps.k 1]),'bits_imag',randi([0 1], [ps.k 1]));
%% decoding
b_rx        = ps.decoding(s);
%% compute the bit erro probability
pe_b        = sum(not(b==b_rx))/length(b);