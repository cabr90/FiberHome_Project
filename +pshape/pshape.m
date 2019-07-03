classdef pshape < handle
    
    
    properties (SetAccess='immutable')
        M       % modulation order
        dm_rate % rate of the dm
        rate    % rate of ps qam
        N       % dm output size
        k       % dm input size
        dm_type % the type of the dm
    end
    
    properties (Access='private')
        dm      % distribution matcher for the ps encoding
    end
    
    methods (Access='public')
        
        function ps               = pshape(dm,M,N,r,varargin)
            % PSHAPE is the class constructor where :
            % -dm : is the distribution matcher type (only ccdm up to now)
            % -M  : is the qam modulation order
            % -N  : is the size of the ccdm output
            % -r  : is the rate of the system
            % -opt: numbits for ess
            
            import pshape.*;
            import pshape.dm.*;
            
            if(strcmp(dm,'ccdm'))
                ps.M       = M;
                ps.dm_rate = (r-2)/2;
                ps.N       = N;
                ps.rate    = r;
                M_DM       = sqrt(M)/2;
                dm         = ccdm(M_DM, N, ps.dm_rate);
                dm.k       = floor(dm.k);
                ps.dm      = dm;
                ps.k       = dm.k;
                ps.dm_type = 'ccdm';
            elseif(strcmp(dm,'ess'))
                ps.M       = M;
                ps.dm_rate = (r-2)/2;
                ps.N       = N;
                ps.rate    = r;
                M_DM       = sqrt(M)/2;
                if isempty(varargin)
                    dm     = ess(M_DM,N,ps.dm_rate);
                elseif strcmp(varargin{1},'numbits')
                    dm     = ess(M_DM,N,rate_DM,varargin{:});
                else
                    error('the name of the extra input must be: numbits');
                end
                dm.k       = round(ps.dm_rate*ps.N);
                ps.dm      = dm;
                ps.k       = dm.k;
                ps.dm_type = 'ess';
            end
            
        end
        
        function [symb,bits,pqam] = encoding(ps,varargin)
            % ENCODING generates a qam sequence of length N encoded with
            % the dm initialized at the constructor of the class
            norm = 1/sqrt(2/3*(ps.M-1));
            if not(isempty(varargin)) && (length(varargin) == 2 || length(varargin) == 6) && strcmp(varargin{1},'UnitAveragePower') && not(varargin{2})
                norm = norm * varargin{2}+1;
            end
            
            
            if isempty(varargin) || length(varargin) == 2 
                % Generate bits
                bits_real      = randi([0 1], [ps.dm.k 1]);     % bits for DM for the real part of the amplitude of the MQAM symbol
                bits_imag      = randi([0 1], [ps.dm.k 1]);     % bits for DM for the imaginary part of the amplitude of the MQAM symbol
            elseif not(length(varargin) == 4 || length(varargin) == 6)
                error('if you want to pass the bits to encode the right way is: ps.encoding(bits_real,[array_bits],bits_imag,[array_bits]) ')
            elseif length(varargin) == 4 && strcmp(varargin{1},'bits_real') && strcmp(varargin{3},'bits_imag')
                bits_real      = varargin{2};                   % bits for DM for the real part of the amplitude of the MQAM symbol
                bits_imag      = varargin{4};                   % bits for DM for the imaginary part of the amplitude of the MQAM symbol
            elseif length(varargin) == 6 && strcmp(varargin{3},'bits_real') && strcmp(varargin{5},'bits_imag')
                bits_real      = varargin{4};                   % bits for DM for the real part of the amplitude of the MQAM symbol
                bits_imag      = varargin{6};                   % bits for DM for the imaginary part of the amplitude of the MQAM symbol
            else
                error('args are: ps.encoding(bits_real,[array_bits],bits_imag,[array_bits])')
            end
            
            bits_sign_real = randi([0 1], [ps.N 1]);        % bits for the sign of the real part of the MQAM symbol
            bits_sign_imag = randi([0 1], [ps.N 1]);        % bits for the sign of the imaginary part of the MQAM symbol
            
            % Generate amplitudes with dm
            if strcmp(ps.dm_type,'ccdm')
                amplitudes_real = 2*ps.dm.encode(bits_real).'-1; % amplitudes for real part of the MQAM
                amplitudes_imag = 2*ps.dm.encode(bits_imag).'-1; % amplitudes for imag part of the MQAM
            else
                amplitudes_real = ps.dm.encode(bits_real).'; % amplitudes for real part of the MQAM
                amplitudes_imag = ps.dm.encode(bits_imag).'; % amplitudes for imag part of the MQAM
            end
            % Generate bit and symbols
            bits     = [bits_sign_real;bits_sign_imag;bits_real;bits_imag];
            symb     = norm*(ps.beta(bits_sign_real).*(amplitudes_real)+1j*ps.beta(bits_sign_imag).*(amplitudes_imag));
            pqam     = ps.get_prob(symb);
            
        end
        
        function [bits]           = decoding(ps,symbols,varargin)
            % DECODING decode the sequence in SYMBOLS using the dm
            % initialized at the class constructor
            norm        = 1/sqrt(2/3*(ps.M-1));
            if not(isempty(varargin)) && length(varargin) == 2 && strcmp(varargin{1},'UnitAveragePower') && not(varargin{2})
                norm = norm * varargin{2}+1;
            end
            
            symbols     = symbols./norm;
            
            symbol_real = real(symbols);
            symbol_imag = imag(symbols);
            
            bits_sgn_real = ps.betainv(symbol_real);
            bits_sgn_imag = ps.betainv(symbol_imag);
            
            if strcmp(ps.dm_type,'ccdm')
                bits_real = ps.dm.decode(round((abs(symbol_real)+1)/2,13));
                bits_imag = ps.dm.decode(round((abs(symbol_imag)+1)/2,13));
            else
                bits_real = ps.dm.decode(round(abs(symbol_real),13));
                bits_imag = ps.dm.decode(round(abs(symbol_imag),13));
            end
            bits=[bits_sgn_real;bits_sgn_imag;bits_real.';bits_imag.'];
            
        end
        
        function [amps,bits]      = amps_encoding(ps,varargin)
            
            if isempty(varargin)
                % Generate bits
                bits      = randi([0 1], [ps.dm.k 1]); % bits for DM for the real part of the amplitude of the MQAM symbol
            elseif not(length(varargin) == 1)
                error('if you want to pass the bits to encode the right way is: ps.encoding("bits",[array_bits]')
            elseif strcmp(varargin{1},'bits')
                bits      = varargin{2};               % bits for DM for the real part of the amplitude of the MQAM symbol
            else
                error('args are: ps.encoding("bits",[array_bits])')
            end
            
            % Generate amplitudes with dm
            if strcmp(ps.dm_type,'ccdm')
                amps = 2*ps.dm.encode(bits).'-1; % amplitudes for real part of the MQAM
            else
                amps = ps.dm.encode(bits).';     % amplitudes for real part of the MQAM
            end
            
        end
        
        function [bits]           = amps_decoding(ps,amps)
            if strcmp(ps.dm_type,'ccdm')
                bits = ps.dm.decode(round((abs(amps)+1)/2,13));
            else
                bits = ps.dm.decode(round(abs(amps),13));
            end
        end
        
    end
    
    methods (Access='private')
        
        function[p]=get_prob(ps,simbols)
            norm       = 1/sqrt(2/3*(ps.M-1));
            m          = ps.M;
            X          = real(pammod((0:sqrt(m)-1),sqrt(m)));
            Y          = real(pammod((0:sqrt(m)-1),sqrt(m)));
            MQAM_stars = norm*(X+1i*Y.');
            p          = arrayfun(@(x)length(find(simbols==x)),MQAM_stars)/length(simbols);            
        end
        
    end
    
    methods (Access='private',Static)
        
        function[sgn]=beta(bits)
            sgn=2*bits-1;
        end
        
        function[bit]=betainv(symbols)
            bit=(sign(symbols)+1)/2;
        end
        
    end
    
    
end