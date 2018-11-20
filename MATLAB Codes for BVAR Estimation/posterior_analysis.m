% posterior_analysis.m
function [imp_resp_upper, imp_resp_lower, pbound_resp, imp_resp_diff_qus] = posterior_analysis(imp_upper, imp_lower, qus, alphas)
[nsave, nvar, ihor] = size(imp_upper);
imp_resp_upper = squeeze(quantile(imp_upper,qus));
imp_resp_lower = squeeze(quantile(imp_lower,qus));
imp_resp_diff = imp_upper - imp_lower;
imp_resp_diff_qus = squeeze(quantile(imp_resp_diff,qus));
mean_diff = mean(imp_resp_upper,1) - mean(imp_resp_lower,1);
% var_upper = var(imp_resp_upper,0,1);
% var_lower = var(imp_resp_upper,0,1);
pooled_std = sqrt(var(imp_resp_upper,0,1)/nsave + var(imp_resp_lower,0,1)/nsave);
criticals = norminv(1-alphas/2);
criticals = sort([-criticals 0 criticals]);
pbound_resp = zeros(length(criticals),nvar,ihor);
for i = 1:ihor
    pbound_resp(:,:,i) = criticals' * squeeze(pooled_std(:,:,i)) + mean_diff(:,:,i);
end
end