function cor_table = correlate_brain_behavior(d, lo, Table, correlation)

measure = {'pain_avg','spon_pain_rating_mean'};
area    = {'m_dorsal_insula','m_posterior_insula','A1','mPFC','precuneus'};
fields  = {'all_s_l','all_s_h'};           % low/high
rowNames = {'clbp_s_l','clbp_s_h'};        % to label rows in the table

results = struct(); 

for n = 1:numel(area)
    aName = area{n};

    for m = 1:numel(measure)
        measName = measure{m};

        % Preallocate as 2x1 (for low/high)
        p_value = nan(numel(fields),1);
        Rho     = nan(numel(fields),1);

        for i = 1:numel(fields)
            fi = fields{i};

            % X: concatenate across patient groups G1..G3
            X = rmoutliers([d.G1.S1.(measName); ...
                 d.G2.S1.(measName); ...
                 d.G3.S1.(measName)]);

            % Y: same groups, given ROI/area field and low/high subfield
            Y = rmoutliers([lo.roi.G1.S1.(aName).(fi); ...
                 lo.roi.G2.S1.(aName).(fi); ...
                 lo.roi.G3.S1.(aName).(fi)]);

            % Ensure column vectors
            X = X(:);
            Y = Y(:);

            % Correlation (Pearson), ignore NaNs
            [rho, pval] = corr(X, Y, 'Rows','complete', 'Type','Pearson');

            p_value(i) = pval;
            Rho(i)     = rho;
        end

        % Build a tidy 2-row table
       T = table(p_value, Rho, ...
                  'VariableNames', {'p_value','Rho'}, ...
                  'RowNames', rowNames);

        % Store it under results.corr.<measure>.<area>.clbp
        results.corr.(measName).(aName) = T;
    end
end


end 