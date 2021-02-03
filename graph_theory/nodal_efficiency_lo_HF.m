function [Eff_thr_nodal, E_norm_nodal] = nodal_efficiency_lo_HF(E)

% ######################################################################
%
% Script to produce the:
%       - Normalized nodal efficiency (E_norm_nodal)
%       - Detected hubs - thresholded (Eff_thr_nodal)
%
% The normalized nodal efficiency is obtained by dividing the nodal
% efficiency by the mean of all nodes
%
% The threshold for hub identification is defined as the mean + SD
%
%   Henrique Fernandes 2013
%   based on: Lo et al 2010
%
% ######################################################################

Eff_thr_nodal=zeros(size(E.mean));
Eff_avg=mean(E.mean);

E_norm_nodal.mean=E.mean/Eff_avg;
E_norm_nodal.stdev=E.stdev/Eff_avg;


Eff_avg_norm=mean(E_norm_nodal.mean);
Eff_norm_stdev=std2(E_norm_nodal.mean);
Eff_thresh=Eff_avg_norm+Eff_norm_stdev;


Eff_nodes_arr=find(E_norm_nodal.mean>Eff_thresh);

for i=1:size(E.mean)
   if ismember(i,Eff_nodes_arr)
%        Eff_thr_nodal(i)=E.mean(i);
       Eff_thr_nodal(i)=E_norm_nodal.mean(i);
   end
end

end




