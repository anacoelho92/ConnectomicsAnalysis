function [pop_ref, Ci_ref, Qa_mean] = get_reference_scheme(mats_pop1, mats_pop2, pop1_label, pop2_label, Ci1, Ci2)

% Function to get reference scheme to compare fingerprints of modular
% connectivity between groups
% Compares the mean score of community structure goodness-of-fit and the
% population with higher score is the reference
%
%
% Inputs:   mats_pop1/pop2: 3D weighted connectivity matrices for each group
%
%           pop1/pop2_label: label of each group
%
%           Ci1/Ci2: community structure of each group
%
% Outputs:  pop_ref: name of reference group
%
%           Ci_ref: reference scheme
%
%
% Ana Coelho 15-04-2020
%

% compute goodness-of-fit of each group
n_pop1 = 1:size(mats_pop1,3);
for p=1:length(n_pop1)
    Qa1(p)=modularity_fitting_und_sign(mats_pop1(:,:,p),Ci1);
end

n_pop2 = 1:size(mats_pop2,3);
for p=1:length(n_pop2)
    Qa2(p)=modularity_fitting_und_sign(mats_pop2(:,:,p),Ci2);
end

% check which group has higher mean score
if mean(Qa1) > mean(Qa2)
    pop_ref = pop1_label;
    Ci_ref = Ci1;
    Qa_mean = mean(Qa1);
else
    pop_ref = pop2_label;
    Ci_ref = Ci2;
    Qa_mean = mean(Qa2);
end


end
