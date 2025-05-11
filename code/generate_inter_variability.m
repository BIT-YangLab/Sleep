function inter_varia = generate_inter_variability(profile_list, interrupt_ths)
% reference: Mueller et al., 2013
    inter_varia = zeros(size(profile_list{1}, 1), 1);
    nc1 = 0;
    for si = 1:length(profile_list)-1
        disp(['compute inter variability: ' num2str(si) '/' num2str(length(profile_list)-1) ])
        for sj = si+1:length(profile_list)
            corr_tmp = paircorr_mod(profile_list{si}', profile_list{sj}');
            nc1 = nc1 + 1;
            inter_varia(:, nc1) = diag(corr_tmp);
        end
        if si > interrupt_ths
            break
        end
    end
    inter_varia = nanmean(1-inter_varia, 2);


end