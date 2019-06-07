classdef ccdm < handle
    
    properties
        k
        n
        M
        pA
        type_seq
        kmax
    end
    
    methods (Access='public')
        
        function obj = ccdm(M, n, rate, varargin)            
            obj.n = n;
            obj.M = M;
            [obj.pA, obj.k, obj.type_seq, obj.kmax] = obj.get_ccdm_params(struct('M', obj.M, 'n', obj.n, 'k', n * rate));
            
        end
              
        function symbols = encode(obj, bits)
            if length(bits) ~= obj.k
                error('Parameter "bits" must have size: %d. Given was: %d', obj.k, length(bits));
            end
            
            counts     = uint64(obj.pA'.*obj.n);
            cum_counts = [0, cumsum(counts)];
            p_src      = [0.5,0.5];
            p_code     = double((cum_counts(2:end)-cum_counts(1:end-1)))/double(cum_counts(end));
            
            Isrc       = [0.0,1.0];
            
            Cand_int   = obj.refine([0.0,1.0],p_code);
            
            len        = obj.n;
            symbols    = zeros(1,obj.n);          
            
            j = 0;
            
            for h = 1:length(bits)
                Isrc               = obj.readSymb(Isrc,p_src,bits(h)+1);
                index              = obj.findCandidate(Isrc,Cand_int);
                while index
                    j = j+1;                 
                    [symbols(j),p_code,len,counts]  = obj.update_code(len,counts,index);                    
                    Icode              = Cand_int(index,:);
                    [Isrc,Icode]       = obj.rescale(Isrc,Icode);
                    Cand_int           = obj.refine(Icode,p_code);
                    index              = obj.findCandidate(Isrc,Cand_int);
                end
            end
            
            if(not(len == 0))
                index   = obj.findCandidateOutput(Isrc,Cand_int,1);
                while (index)
                    if index
                        j = j+1;
                    end
                    [symbols(j),p_code,len,counts]  = obj.update_code(len,counts,index);
                    if len == 0
                        break;
                    end
                    Icode              = Cand_int(index,:);
                    Cand_int           = obj.refine(Icode,p_code);
                    index              = obj.findCandidateOutput(Isrc,Cand_int,1);
                    
                end
            end
            
            
        end
        
        function bits = decode(obj, symbols)
            
            if length(symbols) ~= obj.n
                error('Parameter "symbols" to decode() must have size: %d. Given was: %d', obj.n, length(symbols));
            end
                      
            counts     = uint64(obj.pA'.*obj.n);
            cum_counts = [0, cumsum(counts)];
            p_src      = [0.5,0.5];
            p_code     = double((cum_counts(2:end)-cum_counts(1:end-1)))/double(cum_counts(end));
            
            Cand_int   = obj.refine([0.0,1.0],p_src);
            
            len        = obj.k;
            i          = 1;
            j          = 0;
            bits       = zeros(1,obj.k);
            
            while i<length(symbols)
                h            = i+1;
                [~,p_code_future,~,counts_future] = obj.update_code(len,counts,symbols(i));
                Icode        = obj.readSymb([0,1],p_code,symbols(i));
                scaling      = 0;
                while(~scaling)
                    
                    while obj.findCandidate(Icode,Cand_int) 
                        j           = j + 1;
                        bits(j)     = obj.findCandidate(Icode,Cand_int)-1;
                        Isrc        = Cand_int(bits(j)+1,:);
            
                         
                        if(j >= len)
                            return;
                        end
                        
                        [forward_steps,Isrc,p_code,counts] = obj.look_for_scaling(Isrc,p_code,counts);
                        if forward_steps
                            Cand_int  = obj.refine(Isrc, p_src);
                            i         = i+forward_steps;
                            scaling   = 1;
                            break;
                        else
                            Cand_int    = obj.refine(Isrc, p_src);
                        end
                        
                    end
                    
                    if (h == obj.n)
                        Icode(2) = Icode(1) + 0.01*(Icode(2)-Icode(1));
                    else
                        Icode = obj.readSymb(Icode,p_code_future,symbols(h));
                        [~,p_code_future,~,counts_future] = obj.update_code(len,counts_future,symbols(h));
                        h = h+1;
                    end
                    
                end
            end
            
        end
        
    end
    
    methods (Access='private',Static)
        
        function [pA, k, ni, kmax]  = get_ccdm_params(ccdm)
                 
            [ni_uniform, ~] = idquant(ones(1,ccdm.M)/ccdm.M, ccdm.n);
            kmax            = floor(n_choose_ks_recursive_log2(ccdm.n, ni_uniform));
            
            if ccdm.k > kmax
                error('For this blocklength and modulation, the rate is to high')
            end
            
            ccdm.p_target = exp(-((1:2:2*ccdm.M-1)/ccdm.M).^2);
            ccdm.p_target = ccdm.p_target/sum(ccdm.p_target);
            
            [nu,~,exitflag] = fzero(@(nuVal) find_rate(nuVal, ccdm.p_target, ccdm.n) - ccdm.k/ccdm.n, [1e-6 1e3]);
            
            addidionalbit = 0;
            while(find_rate(nu, ccdm.p_target, ccdm.n) - ccdm.k/ccdm.n) < 0
                addidionalbit = addidionalbit+1;
                nu = fzero(@(nu_var) find_rate(nu_var, ccdm.p_target, ccdm.n) - (ccdm.k+addidionalbit)/ccdm.n, [1e-6 1e3]);
            end
            
            if exitflag ~=1
                error('No distribution scaling can be found to support the desired rate.');
            else
                ccdm.p_target = ccdm.p_target.^nu/sum(ccdm.p_target.^nu);
                [ni, pA] = idquant(ccdm.p_target, ccdm.n);
                k = ccdm.k;
                entropyrate = H(pA);
                rateloss = entropyrate - ccdm.k/ccdm.n;
                if rateloss > 0.1
                    warning('Your current selection will lead to a rateloss greater than 0.1 bit/symbol. Please increase the DM output length to avoid that loss.');
                end
            end
        
            function rate = find_rate(nu, p_target, n)
                
                [n_i, ~] = idquant(p_target.^nu/sum(p_target.^nu), n);
                
                kcalc = floor(n_choose_ks_recursive_log2(n, n_i));
                
                rate = kcalc / n;
                
            end
            
            function entropy = H(pX)
                
                if numel(pX)==1
                    pX = [pX 1-pX];
                end
                
                entropy = -pX .* log2(pX);
                entropy(~isfinite(entropy)) = 0;
                entropy = sum(entropy);
                
            end
            
            function [n_i,p_quant] = idquant(p,n)
                m = length(p);
                n_i = zeros(m,1);
                t = log(1./p);
                p_quant = t;
                
                
                for h=1:n
                    [~,index] = min(p_quant);
                    cj = n_i(index)+1;
                    n_i(index) = cj;
                    p_quant(index) = (cj+1)*log(cj+1)-cj*log(cj)+t(index);
                end
                p_quant = n_i/n;
            end
            
            function [ out ] = n_choose_k_iter_log( n,k )
                if k > n || k < 0
                    error('k must be smaller than n and bigger and non-negative')
                end
                k = min(k,n-k);
                i = 1:k;
                out = sum(log2((n-(k-i))./(i)));
                
            end
            
            function [ out ] = n_choose_ks_recursive_log2( n,k )
                
                out = 0;
                k = sort(k);
                
                for i = 1:length(k)-1
                    out = out + n_choose_k_iter_log(n,k(i));
                    n = n - k(i);
                end
            end
            
        end
        
        function [type_seq] = get_type(symbols)
            symbols_u = unique(symbols);
            type_seq  = zeros(length(symbols_u),1);
            for ii=1:length(symbols_u)
                type_seq(ii) = sum(symbols_u(ii) == symbols);
            end
        end
        
        function [I1,I2] = rescale(I1,I2)
            I1(1) = (I1(1)-I2(1))/[I2(2)-I2(1)];
            I1(2) = (I1(2)-I2(1))/[I2(2)-I2(1)];
            if(I1(2)>1.0)
                I1(2) = 1.0;
            end            
            I2    = [0,1];
        end
        
        function [c,p,l,cn] = update_code(l,counts,s)
            counts(s)  = counts(s)-1;
            cum_counts = [0, cumsum(counts)];
            p     = double((cum_counts(2:end)-cum_counts(1:end-1)))/double(cum_counts(end));
            l     = l-1;
            c     = s;
            cn    = counts;
        end
        
        function idx = findCandidateOutput(Isrc,List,enco)
            if enco
                idx = find(Isrc(1) <= List(:,1) & Isrc(2) > List(:,1) & List(:,1) ~= List(:,2));
                if (isempty(idx))
                    idx = 0;
                end
                idx = idx(1);
            else
                idx = find(List(:,1) <= Isrc(1)  & List(:,2) > (Isrc(2)-1e-6) );
                if (isempty(idx))
                    idx = 0;
                end
                idx = idx(1);
            end
        end
        
        function idx = findCandidate(Isrc,List)
            
            idx = find(List(:,1) <= Isrc(1)  & List(:,2) > (Isrc(2)-1e-13) );
            if (isempty(idx))
                idx = 0;
            end
            idx = idx(1);
            
        end
        
        function Isrc = readSymb(I,p,symb)
            F       = [0,cumsum(p)];
            Isrc(1) = I(1) + (I(2)-I(1))*F(symb);
            Isrc(2) = I(1) + (I(2)-I(1))*F(symb+1);            
        end
        
        function list = refine(I,p)
            F         = [0,cumsum(p)];
            list(1,:) = [I(1) + (I(2)-I(1))*F(1),I(1) + (I(2)-I(1))*F(2)];
            for i = 2:length(p)
                list(i,:) = [I(1) + (I(2)-I(1))*F(i),I(1) + (I(2)-I(1))*F(i+1)];
            end
        end
        
        function   [forward_steps,Isrc,p_code,counts]  = look_for_scaling(Isrc,p_code,counts)
            import pshape.*;
            import pshape.dm.*;
            
            forward_steps = 0;
            ccCandidate = ccdm.refine([0 1], p_code);            
            
            while ccdm.findCandidate(Isrc,ccCandidate)
                index                = ccdm.findCandidate(Isrc,ccCandidate);
                forward_steps        = forward_steps+1;
                [~,p_code,~,counts]  = ccdm.update_code(0,counts,index);
                Icode                = ccCandidate(index,:);
                [Isrc,Icode]         = ccdm.rescale(Isrc,Icode);
                ccCandidate          = ccdm.refine(Icode,p_code);
            end
            
        end
        
    end
    
end
