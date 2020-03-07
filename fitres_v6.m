%% fitres_v6 - batch analysis of radial pulse data using kreservoir_vXX by KHP
%  Copyright 2019 Alun Hughes
%  This software is distributed under under the terms of the GNU General Public License
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% http://www.gnu.org/licenses/gpl.html

%% Versions
% v1 First stable version(08/07/15)
% v2 Suppressed figure display without preventing save, removed circles
% indicating Pn/Tn
% v3 added HRV measures
% v3 replaced fsg1521 with sgolayfilt
% v5 excluded any data at the end of the trace where pressure rises
% v5 included max dp as the cycle marker (in addition to SBP) to test out
% v6 removed HRV to separate function
%%
function  [P_av, Pr_av,Pinf_av,Pn_av,Tn_av,Sn_av,fita_av,fitb_av,rsq_av]=fitres_v6(pulse,sampling_rate) 
%%
% analyse the average beat
% *****NB  peripheral pulse starts at the foot!!!
    p_av=pulse';
    j=find(p_av==0); % gets rid of any zeros which occasionally appear at the end of the trace (probably due to the transfer function-induced shift)   
    p_av=p_av(1:j(1)-1);
  

    % all data
    T_av=length(p_av)/sampling_rate;
 
    % prevent upturn in pressure affecting fit (different methods compared
    % but settled on method a (other methods retained but commented out
    
    % (a) exclude diff (p)>0 in diastole
    diffp=diff(p_av);
    cut=find(diffp<0,1,'last');
    P_av=p_av(1:cut);
    
    % (b) exclude last 40ms (last 6 samples)
    % P_av=p_av(1:end-6);  % in tests identical to (a)
    % P_av=p_av(1:end-18); % rather extreme 
    % (c) include on last 2/3 of diastole [suggested by Leif Rune Hellevik
    % in Cardiovascular Biomechanics] - not done since problem is with end
    % - Pinf not the early part of the fit
    
    % fit average beat using kreservoir
    [Pr_av,fita_av,fitb_av,Pinf_av,Tn_av,Pn_av]=kreservoir_v14(P_av,T_av, sampling_rate);

    % Calculate R^2 (cofficient determination) for Pr fit in diastole
    Sn_av=round(Tn_av*sampling_rate);         % parameters for length diastole for R2
    C = corrcoef(P_av(Sn_av:end),Pr_av(Sn_av:end)); % calculate correlation for diastolic fit
    rsq_av = C(1,2)^2;

  end