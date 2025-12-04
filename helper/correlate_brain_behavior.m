function results = correlate_brain_behavior(d, lo, correlation)


if strcmpi(correlation, 'brain_clinpain')

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
                X = [d.G1.S1.(measName); ...
                    d.G2.S1.(measName); ...
                    d.G3.S1.(measName)];

                % Y: same groups, given ROI/area field and low/high subfield
                Y = [lo.roi.G1.S1.(aName).(fi); ...
                    lo.roi.G2.S1.(aName).(fi); ...
                    lo.roi.G3.S1.(aName).(fi)];

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


elseif strcmpi(correlation, 'brain_unpleasanteness')

    
    measure    = {'acute_mean_sound_lo','acute_mean_sound_hi'};
    area    = {'m_dorsal_insula','m_posterior_insula','A1','mPFC','precuneus'};
    rowNames_clbp = {'clbp_s_l','clbp_s_h'};
    rowNames_hc   = {'hc_s_l','hc_s_h'};
    fields  = {'all_s_l','all_s_h'};           % low/high
    rowNames = {'clbp_s_l','clbp_s_h'};        % to label rows in the table

    clear p_value;clear pval;clear rho;clear Rho;clear Rho_hc; clear p_value_hc;
    area = {'dorsal_insula' 'posterior_insula' 'A1' 'mPFC', 'precuneus'};
    Condition = ["Sound_low";"Sound_high"];
    stimtype = {'acute_mean_sound_lo' 'acute_mean_sound_hi'};

    for n = 1: numel(area)
        %figure("Name",(area{n}))
        for i = 1:2
            j = i + 2;
            X = d.CLBP_metadata.(stimtype{i});
            Y = roi.(area{n}).(fields{j});
            [rho,pval] = corr(X,Y,'Rows','complete', 'type','Pearson');
            p_value(i,:) = pval;
            Rho(i,:) = rho;
            
            p = j + 6;
            X_hc = d.HC_metadata.(stimtype{i});
            Y_hc = roi.(area{n}).(fields{p});
            [rho,pval] = corr(X_hc,Y_hc,'Rows','complete', 'type','Pearson');
            p_value_hc(i,:) = pval;
            Rho_hc(i,:) = rho;

        end
         
            table.corr.stim_pain.(area{n}).clbp =  table(Condition, p_value, adj_p, Rho);
            table.corr.stim_pain.(area{n}).hc   =  table(Condition, p_value_hc,adj_p_hc, Rho_hc);

    end
    results = struct();

    for n = 1:numel(area)
        aName = area{n};  %#ok<NASGU> % kept for symmetry; not used below but available if you later need it

        for m = 1:numel(measure)
            measName = measure{m};

            % ---- CLBP (G1..G3 pooled) ----
            p_value_clbp = nan(numel(unpleasant),1);
            rho_clbp     = nan(numel(unpleasant),1);

            % common X for CLBP: concatenate across patient groups
            X_clbp = [d.G1.S1.(measName); d.G2.S1.(measName); d.G3.S1.(measName)];
            X_clbp = X_clbp(:);

            for i = 1:numel(unpleasant)
                ui = unpleasant{i};

                % Y for CLBP: same order/subjects as X_clbp
                Y_clbp = [d.G1.S1.(ui); d.G2.S1.(ui); d.G3.S1.(ui)];
                Y_clbp = Y_clbp(:);

                % pairwise-complete mask (no outlier desync!)
                keep = isfinite(X_clbp) & isfinite(Y_clbp);
                if nnz(keep) < 3
                    rho = NaN; pval = NaN;
                else
                    [rho, pval] = corr(X_clbp(keep), Y_clbp(keep), 'Type','Pearson');
                end

                rho_clbp(i)     = rho;
                p_value_clbp(i) = pval;
            end

            T_clbp = table(p_value_clbp, rho_clbp, ...
                           'VariableNames', {'p_value','Rho'});
            T_clbp.Properties.RowNames = rowNames_clbp;

            % ---- HC only ----
            p_value_hc = nan(numel(unpleasant),1);
            rho_hc     = nan(numel(unpleasant),1);

            X_hc = d.HC.S1.(measName);
            X_hc = X_hc(:);

            for i = 1:numel(unpleasant)
                ui = unpleasant{i};
                Y_hc = d.HC.S1.(ui);
                Y_hc = Y_hc(:);

                keep = isfinite(X_hc) & isfinite(Y_hc);
                if nnz(keep) < 3
                    rho = NaN; pval = NaN;
                else
                    [rho, pval] = corr(X_hc(keep), Y_hc(keep), 'Type','Pearson');
                end

                rho_hc(i)     = rho;
                p_value_hc(i) = pval;
            end

            T_hc = table(p_value_hc, rho_hc, ...
                         'VariableNames', {'p_value','Rho'});
            T_hc.Properties.RowNames = rowNames_hc;

            % Store
            results.behavior.(measName).clbp = T_clbp;
            results.behavior.(measName).hc   = T_hc;
        end
    end
end

end


