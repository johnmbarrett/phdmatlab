function [I,Icont,Icount,Ipart] = binlessInfo(s,R,varargin)
    options = getopt('kld=false base=2 starttime=0 triallength=1.0 label=''ch1'' samplingfrequency=25000 stratification_strategy=2 singleton_strategy=0 min_embed_dim=1 max_embed_dim=2 warpfun=NaN embedfun=NaN',varargin{:});
    
    assert(size(R,1) == numel(s),'R must have as many rows as there are trials');
    
    if ~iscell(R)
        sizeR = num2cell(size(R));
        R = mat2cell(R,ones(sizeR{1},1),sizeR{2:end});
    end
    
%     S = makeSTAToolkitStruct(s,R,options.triallength,options.label,options.samplingfrequency,options.starttime);
    
    warningState = warning('query');
    warning off; %#ok<WNOFF>
    
%     [data,counts,categories] = binlessopen(S);
    
    data = R;
    counts = int32(cellfun(@numel,data));
    [uniqueCategories,~,categories] = unique(s);
    categories = int32(categories-1);
    
%     warped = binlesswarp(data);

    m = size(data,2);
    warped = cell(size(data));
    
    for hh = 1:m
        isRow = cellfun(@isrow,data(:,hh));
        isCol = cellfun(@iscolumn,data(:,hh));
        cumCounts = [0; cumsum(double(counts(:,hh)))];

        if all(isRow)
            allTimes = [data{:,hh}];
        elseif all(isCol)
            allTimes = vertcat(data{:,hh})';
        else
            allTimes = zeros(1,sum(counts(:,hh)));

            for ii = 1:size(data,1)
                allTimes(cumCounts(ii)+1:cumCounts(ii+1)) = data{ii,hh};
            end
        end 

        if isa(options.warpfun,'function_handle')
            warpFun = options.warpfun;
        else
            warpFun = @(ts) 2*(tiedrank(ts)-0.5)/numel(ts)-1;
        end

        allTimes = warpFun(allTimes);

        for ii = 1:size(data,1)
            warped{ii,hh} = allTimes(cumCounts(ii)+1:cumCounts(ii+1));
        end
    end
    
%     embedded = binlessembed(warped,struct('max_embed_dim',options.max_embed_dim));
    
    d = options.max_embed_dim;
    
    embedded = zeros(size(warped,1),d+1,m);
    
    % TODO : this has stopped working, investigate why
%     for hh = 1:m
%         embedded(:,:,hh) = binlessembed(warped(:,hh),struct('max_embed_dim',d));
%     end
    
    if isa(options.embedfun,'function_handle')
        embedFun = options.embedfun;
    else
        P = generateLegendrePolynomials(d);
        embedFun = @(h,tks) sqrt(2*h+1)*sum(polyval(P{h+1},tks));
    end
    
    for hh = 1:m
        for ii = 1:size(warped,1);
            for jj = 0:d
                if isempty(warped{ii})
                    continue;
                end

                embedded(ii,jj+1,hh) = embedFun(jj,warped{ii,hh});
            end
        end
    end
    
    ns = numel(uniqueCategories);
    
    % TODO : this basically cut & paste from binlessinfom
    if all(logical(options.kld(:)))
        % TODO : m > 1
        assert(m == 1,'Binless KLD is only supported for one-dimensional data at present');
        I = zeros(ns,2);
        
        strata = options.stratification_strategy;
        
        if numel(strata) == 1
            switch strata
                case 0
                    strata = ones(size(counts));
                case 1
                    strata = counts;
                case 2
                    strata = min(counts,d);
            end
        end
        
        [~,iCounts,iStrata] = unique(strata);
        nStrata = max(iStrata);
        
        degenerate = (counts(iCounts) == 0)';
        
        parity = mod((1:size(embedded,1))',2);
        
        Q = cell(1,nStrata);
        
        for hh = 1:2
            for ii = 1:nStrata
                Q{ii} = embedded(iStrata == ii & parity == hh-1,2:end);
            end

            for ii = 1:ns
                P = cell(1,nStrata);

                for jj = 1:nStrata
                    P{jj} = embedded(categories == ii-1 & iStrata == jj & parity == 2-hh,2:end);
                end
                
                I(ii,hh) = mixedKLDivergence(P,Q,degenerate,1,'both',false);
            end
        end
        
        I = max(0,mean(I,2));
        Icont = NaN;
        Icount = NaN;
        Ipart = NaN;
        return;
    end
    
%     [Ipart,Icont,Icount,I] = binlessinfo(embedded,counts,categories,S.M,struct(varargin{:}));
    [Ipart,Icont,Icount,I] = binlessinfom(embedded,counts,categories,ns,varargin{:});
    
%     assert(max(abs([Ipart.value - Ipart Icont - Icont Icount.value - Icount I.value - I])) < 10e-10);
    
    b = log(2)/log(options.base);
    
    Ipart = b*Ipart; %.value;
    Icont = b*Icont;
    Icount = b*Icount; %.value;
    I = b*I; %.value;
    
    warning(warningState);
end
