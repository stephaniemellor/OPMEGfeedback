function [cD] = merge_feedback_Dobjects(D1,D2)
% 

if isfile(fullfile(D1.path, ['c', D1.fname]))
    cD = spm_eeg_load(fullfile(D1.path, ['c', D1.fname]));
else
    S = [];
    S.D = {D1, D2};
    cD = spm_eeg_merge(S);
    save(cD);
end

end