function print_effectsizes(T, label)
% printEffectSizes - Compact effect size printing for LMM results
%
%   T     = table with variables: value, group, intensity, measure
%   label = string/char, e.g. 'ROI' or 'MVPA'

    allMeasures = cellstr(unique(T.measure));
    
    fprintf('\n==== EFFECT SIZES (%s) ====\n', upper(label));
    fprintf('%-16s | %-35s | %-35s\n', label, 'Low Intensity', 'High Intensity');
    fprintf('%s\n', repmat('-',1,90));
    
    for i = 1:numel(allMeasures)
        M  = allMeasures{i};
        Tk = T(T.measure==M,:);
        if isempty(Tk), continue; end

        % Low intensity
        xL = Tk.value(Tk.group=="HC"  & Tk.intensity=="Low");
        yL = Tk.value(Tk.group=="CBP" & Tk.intensity=="Low");
        efL = mes(xL,yL,'hedgesg');
        gL  = efL.hedgesg;  ciL = efL.hedgesgCi(:)';
        mHC_L = mean(xL,'omitnan'); sHC_L = std(xL,'omitnan');
        mCBP_L = mean(yL,'omitnan'); sCBP_L = std(yL,'omitnan');

        % High intensity
        xH = Tk.value(Tk.group=="HC"  & Tk.intensity=="High");
        yH = Tk.value(Tk.group=="CBP" & Tk.intensity=="High");
        efH = mes(xH,yH,'hedgesg');
        gH  = efH.hedgesg;  ciH = efH.hedgesgCi(:)';
        mHC_H = mean(xH,'omitnan'); sHC_H = std(xH,'omitnan');
        mCBP_H = mean(yH,'omitnan'); sCBP_H = std(yH,'omitnan');

        % Print compact line
        fprintf(['%-16s | HC %.2f±%.2f, CBP %.2f±%.2f, g=%+.2f [%+.2f, %+.2f] ', ...
                 '| HC %.2f±%.2f, CBP %.2f±%.2f, g=%+.2f [%+.2f, %+.2f]\n'], ...
            M, mHC_L,sHC_L,mCBP_L,sCBP_L,gL,ciL(1),ciL(2), ...
               mHC_H,sHC_H,mCBP_H,sCBP_H,gH,ciH(1),ciH(2));
    end
end

