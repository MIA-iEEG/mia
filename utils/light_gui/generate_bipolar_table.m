function m_table_effect = generate_bipolar_table(m_table_all, inclusive)

    tmp_table = {};
    subj_elec = cellfun(@(s,e) [s '_' e], m_table_all(:,1), m_table_all(:,2), 'UniformOutput', false);
    [unique_pairs, ~, group_id] = unique(subj_elec);

    for g = 1:numel(unique_pairs)
        sub_table = m_table_all(group_id == g, :);

        % Trier les lignes par index croissant
        indices = cellfun(@double, sub_table(:,3));
        [~, sort_idx] = sort(indices);
        sub_table = sub_table(sort_idx, :);

        for i = 2:size(sub_table, 1)
            row1 = sub_table(i-1, :);
            row2 = sub_table(i, :);

            region1 = row1{5};
            region2 = row2{5};

            montage_name = sprintf('%s_%s%d-%s%d', ...
                row1{1}, row1{2}, row1{3}, row2{2}, row2{3});

            if ~inclusive
                if isequal(region1, region2)
                    new_row = [row1(1:4), {region1}, {montage_name}];
                    tmp_table(end+1, :) = new_row;
                else
                    new_row1 = [row1(1:4), {region1}, {montage_name}];
                    new_row2 = [row2(1:4), {region2}, {montage_name}];
                    tmp_table(end+1, :) = new_row1;
                    tmp_table(end+1, :) = new_row2;
                end
            else
                % inclusive: inclure toutes régions distinctes
                regions = unique({region1, region2});
                for r = 1:length(regions)
                    region_label = regions{r};
                    new_row = [row1(1:4), {region_label}, {montage_name}];
                    tmp_table(end+1, :) = new_row;
                end
            end
        end
    end

    % Supprimer doublons : créer une clé par ligne (chaîne concaténée)
    row_keys = cellfun(@(row) strjoin(cellfun(@(x) char(string(x)), row, 'UniformOutput', false), '|'), ...
                       num2cell(tmp_table, 2), 'UniformOutput', false);
    [~, ia] = unique(row_keys, 'stable');
    m_table_effect = tmp_table(ia, :);

end
