%% hrv_v1 - batch analysis of radial pulse data using kreservoir_vXX by KHP
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
% v1 First stable version(11/01/20)
% v1.1 Fixed Warning: Colon operands must be real scalars. This warning will become an error
% in a future release. 
%%
function  [P_all,sdsbp_mmhg, nbeats, av_rr_ms, av_rrS_ms, sdnn_ms,...
    sdnnS_ms, rmssd_ms, rmssdS_ms, brs_ms_mmhg,sysloc, dialoc]=hrv_v1(pulse,signal,sampling_rate) 

% analyse the whole signal
    p=signal';                  % Non-calibrated pressure trace (calibrate later)
%     j=find(p==0);               % gets rid of any zeros which occasionally appear at the end of the trace (probably due to the transfer function-induced shift)   
%     p=p(1:j(1)-1);
    
    
    p_av=pulse';                % calibrated pressure trace
    % L=length(p);
    pnorm=(p-min(p))/(max(p)-min(p));
    t=(0:length(p)-1)/sampling_rate;    
    
%     % find systolic peaks
%     [syspeaks, sysloc]=findpeaks(pnorm,'MinPeakDistance',64); % enforces that HR < 120bpm
%     % find diastolic peaks
%     [~, dialoc]=findpeaks(abs(pnorm-1),'MinPeakDistance',64);
%     diapeaks=pnorm(dialoc);
    promfact=0.08; % was 0.25
    [syspeaks, sysloc]=findpeaks(p, 'MinPeakProminence',max(p)*promfact,'MinPeakDistance',50);
    promfact=0.25;
    p_upside=-p-min(-p);
    [~, dialoc] = findpeaks(p_upside, 'MinPeakProminence',max(p_upside)*promfact,'MinPeakDistance',30);
    clear p_upside;
    diapeaks=pnorm(dialoc);

    % calculate SDNN and RMSSD as measures of HRV (sdnn_ms rmssd_ms) - this
    % can be done on the uncalibrated data
    nbeats=length(dialoc);
    % use differences in systolic peaks to determine HRV
    % RR intervals
    rr=diff(sysloc/sampling_rate*1000);    
    % mean RR interval
    av_rr_ms=mean(diff(dialoc/sampling_rate*1000));
    av_rrS_ms=mean(diff(sysloc/sampling_rate*1000));
    % SDnn
    sdnn_ms=std(diff(dialoc/sampling_rate*1000));
    sdnnS_ms=std(diff(sysloc/sampling_rate*1000));
    % RMSSD
    rmssd_ms=rms(diff(diff(dialoc/sampling_rate*1000)));
    rmssdS_ms=rms(diff(diff(sysloc/sampling_rate*1000)));
   
    % calibrate p
    avsbp=max(p_av); avdbp=min(p_av); avpp= avsbp-avdbp;
    avSsig=mean(syspeaks); avDsig=mean(diapeaks); avPPsig=avSsig-avDsig;
    P_all=((pnorm-mean(diapeaks))*avpp)+avdbp;
    sbp=P_all(sysloc);
    dbp=P_all(dialoc);
    sdsbp_mmhg= std(sbp(1:end-1)); % Standard deviation of SBP, mmHg   
  
    % use SBP and rr (next beat) change to estimate BRS 
    diffrr=diff(rr(2:end));             % difference in RR (beat i+1)
%     if length(sbp)!=length(diffrr)
%         % diffrr=diff(rr(2:end-1));
%         diffrr=diff(rr(2:end)); 
%     end    
    diffsbp=diff(sbp(1:end-1)); % difference in SBP (beat i)
    diffsbp(abs(diffsbp)<1) = 0; % to enforce difference of at least 1mmHg
    % find running sequence of 3 increases or decreases in SBP
    signdr=sign(diffrr);
    signds=sign(diffsbp);
    if length(signdr)>length(signds)
        signdr=signdr(1:length(signds));
    elseif length(signds)>length(signdr)
        signds=signds(1:length(signdr));  
    end
    signsp=signdr.*signds;
    seq=[1 1 1];
% Function to calculate baroreflex sensitivity (BRS) in ms/mmHg
if length(signsp) < 3  % Ensure there are at least 3 elements in the sequence
    brs_ms_mmhg = NaN;
else
    % Find indices of positive and negative sequences
    posseq = strfind(signsp, seq);
    negseq = strfind(-signsp, seq);
    
    % Initialize BRS value
    brs_ms_mmhg = NaN;
    
    % Check for positive sequence match
    if ~isempty(posseq)
        % Compute regression coefficient for positive sequence
        indices = posseq(1):(posseq(1) + 2); % Only consider the first match
        coeff = polyfit(abs(diffsbp(indices)), abs(diffrr(indices)), 1);
        brs_ms_mmhg = abs(coeff(1));
    elseif ~isempty(negseq)
        % Compute regression coefficient for negative sequence
        indices = negseq(1):(negseq(1) + 2); % Only consider the first match
        coeff = polyfit(abs(diffsbp(indices)), abs(diffrr(indices)), 1);
        brs_ms_mmhg = abs(coeff(1));
    end
end
