function H = trans_func(filter,f) %#codegen
%% TRANS_FUNC
% This function creates the transfer function of a given filter
% H = trans_func(filter,f)
%
%   *** INPUT VARIABLES ***
%  filter:       filter parameters (shape, order, bandwidth, ... )
%  f:            frequency vector (same normalization as bandwidth)
%

%% Supergaussian filter
%  filter.bw:  low-pass 3dB-bandwidth
%  filter.ord: supergaussian order (1: Gaussian)
%
if strcmp(filter.shape,'G')
    fn=f/filter.bw;         % normalized frequency vector
    H=exp(-0.5*log(2)*((fn).^(2*filter.ord)));

%% Raised-cosine (RC) filter
%  filter.rate:  symbol rate (1/T, for Nyquist condition)
%  filter.roff: roll-off factor
%
elseif strcmp(filter.shape,'RC')
    fn=f/filter.rate;       % normalized frequency vector
    ro=filter.roff;         % roll-off factor
    H=ones(size(fn));
    ind=find((abs(fn)>0.5*(1-ro))&(abs(fn)<0.5*(1+ro)));
    H(ind)=0.5*(1+cos(pi/ro*(abs(fn(ind))-0.5*(1-ro))));
    H(abs(fn)>=0.5*(1+ro))=0;

%% Root-raised-cosine (RRC) filter
%  filter.rate:  symbol rate (1/T, for Nyquist condition)
%  filter.roff: roll-off factor
%
elseif strcmp(filter.shape,'RRC')
    fn=f/filter.rate;       % normalized frequency vector 
    ro=filter.roff;         % roll-off factor
    H=ones(size(fn));
    ind=find((abs(fn)>0.5*(1-ro))&(abs(fn)<0.5*(1+ro)));
    H(ind)=sqrt(0.5*(1+cos(pi/ro*(abs(fn(ind))-0.5*(1-ro)))));
    H(abs(fn)>=0.5*(1+ro))=0;
    
%% Wave-selective switch (WSS) filter, finisar paper
%  filter.bw:  low-pass 3dB-bandwidth
%  filter.ord: supergaussian order (1: Gaussian)
%
elseif strcmp(filter.shape,'WSS')
    fn   = f;       % normalized frequency vector 
    sgr2 = filter.BWotf/(2*sqrt(log(2)));
    B2   = 0.5*filter.B;
    %H=0.5*sgr2*sqrt(pi)*(erf((B2-f)/sgr2)-erf((-B2-f)/sgr2));
    H=0.5*(erf((B2-fn)/sgr2)-erf((-B2-fn)/sgr2));

else
    H=zeros(size(f));
end

