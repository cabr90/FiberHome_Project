classdef ess < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    properties (Access = 'private')
        Emax
        numbits
        S_bits
        NLev
        k_DM
    end
    
    properties (Access = 'public')
        k
        n
        M
        tne
    end
    
    methods (Access='public')
        
        function obj       = ess(M,n,rate,varargin)
            %UNTITLED3 Construct an instance of this class
            %   M is the size of the alphabeth
            %   n is the size of the block
            %   rate is the rate of the dm
            %   optional:
            %      numbits bits per cell in the array
            
            if isempty(varargin)
                numbits = 52;
            elseif strcmp(varargin{1},'numbits')
                numbits = varargin{2};
            else
                error('the name of the extra input must be: numbits');
            end
            
            Emax    = n;
            NLev    = M;
            k_DM    = round(rate.*n);    % number of bits encoded by DM
            S_leq   = 1;                 % dipendera da varargin in futuro quando si introdurra anche la versione ottima
            S_bits  = ceil(k_DM/numbits);
            if ceil((k_DM+1)/numbits)>S_bits
                error('change block length');
            end
            % private properties         
            obj.Emax    = Emax;
            obj.numbits = numbits;
            obj.S_bits  = S_bits;
            obj.NLev    = NLev;
            obj.k_DM    = k_DM;
            
            Nseq        = int64(zeros(1,S_bits)); %2^k
            Nseq(1)     = int64(2^(k_DM-(S_bits-1)*numbits));
            
            [tne,T0e] = obj.Tne(n,Emax,NLev,S_leq,S_bits,numbits);
            while obj.vectint_isless(T0e,Nseq)
                Emax = Emax+8;
                [tne,T0e] = obj.Tne(n,Emax,NLev,S_leq,S_bits,numbits);
            end
            % public properties
            obj.tne = tne;
            obj.k   = k_DM;
            obj.n   = n;
            obj.M   = M;
            
%             save(filename,'tne')
        end
                
        function symbols   = encode(obj,bits)
            
            nbits  = obj.numbits;
            sbits  = obj.S_bits;
            nlev   = obj.NLev;
            
            index    = obj.bits2vectint(reshape(bits,1,[]),nbits,sbits);
            symbols  = obj.Ind2Symb(index,obj.tne,0,sbits,nlev,nbits);
            
        end
        
        function bits      = decode(obj,symbols)
            
            nbits  = obj.numbits;
            sbits  = obj.S_bits;
            nlev   = obj.NLev;
            kdm    = obj.k_DM;
            
            index   = obj.Symb2Ind(symbols,obj.tne,nlev,sbits,nbits);
            bits    = obj.vectint2bits(index,nbits,kdm);
            
        end
        
    end
    
    methods (Access='private',Static)
        
        function [Tne,T0e]  = Tne(N,Emax,NLev,S_leq,S_bits,numbits)
            % input:
            %   N:   block/sequence length
            %   Emax:   maximum energy
            %   NLev: number of levels {1,3,..,2*NLev-1}
            %   S_leq: if = 1 (standard) allows the path with energy <= Emax, if = 0 only those with energy = Emax
            %   S__bits : number of vector components
            
            % The function builds the matrix Tne having n_col rows and N columns
            % the rows represents the possible energy values
            % the columns represent the time/symbols position
            % the matrix Tne contains the number of sequences from the node (n,e) to the final node
            
            import pshape.*;
            import pshape.dm.*;
            
            n_rows = max(1,floor((Emax-N)/8)+1);
            Tne = int64(zeros(n_rows,N,S_bits));
            
            % symbols from 1,3,5,7
            Lev = (1:2:NLev*2-1).';
            
            % Define T_n^e, contains the number of possible sequences from n to N
            % initialization on the last column n=N
            if S_leq == 1
                % all energies less or equal to Emax are considered
                Tne(:,N,S_bits) = int64(ones(n_rows,1));
            else
                % only Emax is taken
                Tne(1,N,S_bits) = int64(ones(1,1));
                %Tne(2:end,N,:) = int64(zeros(n_rows-1,1));
            end
            
            for n = N-1:-1:1
                for k = 1:n_rows
                    for kk = 1:NLev
                        ind = k+ (1-Lev(kk).^2)/8;
                        if ind>=1 && ind<=n_rows
                            a = zeros(numel(ind),S_bits);
                            for kkk=1:numel(ind)
                                a(kkk,:) = Tne(ind(kkk),n+1,:);
                            end
                            Tne(k,n,:) = ess.vectint_sum(Tne(k,n,:),ess.vectint_sumv(a,numbits),numbits);
                        end
                    end
                end
            end
            
            % Define T_0^e, the total number of possible sequences
            T0e = int64(zeros(1,S_bits));
            for kk = 1:NLev
                ind = n_rows+ (1-Lev(kk).^2)/8;
                if ind>=1 && ind<=n_rows
                    a = zeros(numel(ind),S_bits);
                    for kkk=1:numel(ind)
                        a(kkk,:) = Tne(ind(kkk),1,:);
                    end
                    T0e = ess.vectint_sum(T0e,ess.vectint_sumv(a,numbits),numbits);
                end
            end
            
        end
        
        function [symbol]   = Ind2Symb(index,Tne,ind0,S_bits,NLev,numbits)
            % input:
            %   index:
            %   Tne: matrix containing the values of the trellis T_n^e
            %   ind0: starting value for the index (standard = 0)
            %   S_bits: number of vector components
            %   NLev: number of levels {1,3,..,2*NLev-1}
            
            % output:
            %   symobol: symbols sequence
            
            import pshape.*;
            import pshape.dm.*;
            
            % symbols from 1,3,5,7
            Lev2 = ((1:2:NLev*2-1).^2).';
            N   = size(Tne,2);
            n_rows   = size(Tne,1);
            
            symbol2 = zeros(1,N);
            index = ess.vectint_diff(index,ind0,numbits);
            E = 0;
            for n = 1:N
                for m = NLev:-1:1
                    if m>1
                        vect_s = Lev2(1:m-1)+E;
                        vect_ind = n_rows+(n-vect_s)/8;
                        vect_ind_n = vect_ind(vect_ind<=n_rows);
                        vect_ind_n = vect_ind_n(vect_ind_n>=1);
                        a = zeros(numel(vect_ind_n),S_bits,'int64');
                        for kkk=1:numel(vect_ind_n)
                            a(kkk,:) = Tne(vect_ind_n(kkk),n,:);
                        end
                        ind_try = ess.vectint_sumv(a,numbits);
                    else
                        ind_try = zeros(1,S_bits,'int64');
                    end
                    
                    if ~ess.vectint_isless(index,ind_try)
                        symbol2(n) = Lev2(m);
                        index = ess.vectint_diff(index,ind_try,numbits);
                        break
                    end
                end
                E = E + symbol2(n);
            end
            symbol = sqrt(symbol2);
        end
        
        function [index]    = Symb2Ind(symbol,Tne,NLev,S_bits,numbits)
            % input:
            %   symobol: symbols sequence
            %   Tne: matrix containing the values of the trellis T_n^e
            %   NLev: number of levels {1,3,..,2*NLev-1}
            
            % output: index
            
            import pshape.*;
            import pshape.dm.*;
            
            Lev2 = ((1:2:NLev*2-1).^2).';
            N = length(symbol);
            n_rows = size(Tne,1);
            symbol2 = symbol.^2;
            index = zeros(1,S_bits,'int64');
            E = 0;
            
            for n = 1:N
                m = (symbol(n)+1)/2;
                if m>1
                    vect_s = Lev2(1:m-1)+E;
                    vect_ind = n_rows+(n-vect_s)/8;
                    vect_ind_n = vect_ind(vect_ind<=n_rows);
                    vect_ind_n = vect_ind_n(vect_ind_n>=1);
                    a = zeros(numel(vect_ind_n),S_bits,'int64');
                    for kkk=1:numel(vect_ind_n)
                        a(kkk,:) = Tne(round(vect_ind_n(kkk)),n,:);
                        % a(kkk,:) = Tne(vect_ind_n(kkk),n,:);
                    end
                    index = ess.vectint_sum(ess.vectint_sumv(a,numbits),index,numbits);
                end
                E = E + symbol2(n);
            end
            
        end
                
        function [intvec]   = bits2vectint(b,elle,S_bits)
            % generates the decimal number corresponding to the bits b
            k = numel(b);
            intvec = int64(zeros(1,S_bits));
            for kk = S_bits:-1:2
                intvec(kk) = int64(bi2de(fliplr(b(max(1,k-(S_bits-kk+1)*elle+1):numel(b)))));
                b(max(1,k-(S_bits-kk+1)*elle+1):numel(b))=[];
            end
            if numel(b)>1
                intvec(1) = int64(bi2de(fliplr(b)));
            end
        end
        
        function [bits]     = vectint2bits(a,elle,k)
            % generates the bits b corresponding to the decimal number a
            S_bits = numel(a);
            bits = [];
            for kk = S_bits:-1:2
                bcurr = fliplr(de2bi(a(kk)));
                bits = [zeros(1,elle-numel(bcurr)),bcurr,bits];
            end
            bcurr = fliplr(de2bi(a(1)));
            bits = [bcurr,bits];
            bits = double([zeros(1,k-numel(bits)),bits]);
        end
        
        function [diff]     = vectint_diff(a,b,elle)
            % performs the operation a minus b,
            % a and b are row vectors, less than 2^elle
            % d = a-b is row vector, less than 2^elle
            na = numel(a);
            nb = numel(b);
            n = max(na,nb);
            if na<nb
                a = [zeros(1,n-na),a];
            elseif na>nb
                b = [zeros(1,n-nb),b];
            end
            diff = int64(zeros(size(a)));
            r_old = int64(0);
            for kk = n:-1:1
                s_prov = a(kk)-b(kk)-r_old;
                r = int64(abs(floor(double(s_prov)/2^(elle))));
                diff(kk) = s_prov+r*2^elle;
                r_old = r;
            end
        end
        
        function [ind]      = vectint_isless(a,b)
            % ind = 1 if a<b, else ind = 0
            ind = 0;
            if numel(a)<numel(b)
                ind = 1;
            elseif numel(a)==numel(b)
                for k = 1:numel(a)
                    if a(k)<b(k)
                        ind=1;
                        return
                    elseif a(k)>b(k)
                        ind=0;
                        return
                    end
                end
            end
            
        end
        
        function [s]        = vectint_sum(a,b,elle)
            % performs the sum between a and b,
            % a and b are row vectors, less than 2^elle
            % s=a+b is row vector, less than 2^elle
            
            na = numel(a);
            nb = numel(b);
            n = max(na,nb);
            if na<nb
                a = [zeros(1,n-na,'int64'),a];
            elseif nb<na
                b = [zeros(1,n-nb,'int64'),b];
            end
            s = zeros(1,n,'int64');
            r = int64(0);
            for k = n:-1:1
                s_prov = a(k)+b(k)+r;
                r = int64(floor(double(s_prov)/2^(elle)));
                s(k) = s_prov - r*2^(elle);
            end
            if r >0
                s = [r,s];
            end
        end
        
        function [s]        = vectint_sumv(v,elle)
            
            % performs the sum between the rows of v,
            % each component of each row of v is less than 2^elle
            % s=v(1,:)+v(2,:)+... is row vector, less than 2^elle
            
            % input:
            % v is a matrix of nv x n elements
            % each rows contains an int64 number with n components
            % output:
            % s : vector of int64 numbers representing a+b
            n = size(v,2);
            s = int64(zeros(1,n));
            r = int64(0);
            for k = n:-1:1
                s_prov = sum(v(:,k),'native')+r;
                r = int64(floor(double(s_prov)/2^(elle)));
                s(k) = s_prov - r*2^(elle);
            end
            if r>0
                s = [r,s];
            end
end
               
    end
end

