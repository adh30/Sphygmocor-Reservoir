%% kreservoir13  -  calculate reservoir pressure from pressure waveform
%  Copyright 2008 Kim H Parker with some additions from Alun Hughes [ADH] 
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
%
%  http://www.gnu.org/licenses/gpl.html
%% Versions
%  v.03 (06/03/08)
%  v.04 (27/03/08) eliminate anonymous functions for backward compatability
%  v.05 (11/04/08) separate pinf for systole and diastole
%  v.06 (26/04/08) clean up code for use on sphygmocor data sets
%  v.07 (30/07/08) optimise fit for the whole of diastole
%  v.07.1 (05/09/08) fit to model results during diastole rather than data
%  v.07.2 (01/10/08) fixed error in fitting routine
%  v.08 (02/11/09) use polynomial approximation for inverting R(B)
%  v.08 (03/11/09) if convergence fails set pr=min(p)
%  v.09 (18/11/09) use max -dp/dt as marker of diastole rather than notch
%  v.10 (21/01/10) redefine Pres using diastolic exponential decay
%  v.10 (10/02/10) use Simpson's rule to calculate integrals
%  v.11 (31/01/14) various minor corrections to run with cafe_data_v5 [ADH]
%  v.12 (08/07/15) various minor corrections to run with batch data (batch_res_v7 [ADH])
%  v.13 (17/03/19) to run with with batch data (batch_res_v12 [ADH])
%  v.14 (11/01/20) added error trap to calculate b using linear fit to
%  log(P) in diastole if method of moments fails [ADH]
%
%%
function [Pr,A,B,Pinf,Tn,Pn]=kreservoir_v14(P,Tb, sampling_rate)
%   inputs  P    - pressure starting at diastolic pressure
%           Tb   - duration of beat
%           Tn   - time of start of diastole (if missing it is calculated
%                   found by the function; if ==0 it is done interactively
%   outputs Pr   - reservoir pressure waveform
%           A    - rate constant relating Pw to U
%           B    - rate constant = 1/tau
%           Pinf - pressure asymptote
%           Tn   - detected time of the dichrotic notch
%           Pn   - pressure at dichrotic notch

sP=size(P);
if (sP(1)>sP(2))
    P=P';
end
p=P;
Nb=length(P);
dt=Tb/(Nb-1);
t=(0:dt:Tb);

% find the minimum dp to determine Tn
    dp=fsg721(p);
    [~,nn]=min(dp);
    tn=t(nn);
 
% calculate moments of pressure during diastole using model d=a*exp(-bt)+c
pd=p(nn:end);       % pressure during diastole
td=t(nn:end)-tn;    % number of samples in diastole
Td=Tb-tn;           % duration of diastole
dt=td(2)-td(1); 
N=length(pd)-1;     % length of diastole-1

% calculation of E2/E1 using Simpson's rule
E0=(pd(1)+4*sum(pd(2:2:N))+2*sum(pd(3:2:N))+pd(N+1))*dt/(3*Td);
if (mod(N,2))
    E0=E0+(pd(N)+pd(N+1))*dt/(6*Td);
end
pde1=(pd-E0).*exp(td/Td);
E1=(pde1(1)+4*sum(pde1(2:2:N))+2*sum(pde1(3:2:N))+pde1(N+1))*dt/3;
if (mod(N,2))
    E1=E1+(pde1(N)+pde1(N+1))*dt/6;
end
pde2=(pd-E0).*exp(2*td/Td);
E2=(pde2(1)+4*sum(pde2(2:2:N))+2*sum(pde2(3:2:N))+pde2(N+1))*dt/3;
if (mod(N,2))
    E2=E2+(pde2(N)+pde2(N+1))*dt/6;
end
r=E2/E1;

% invert R(BTd)=r
global Rexp_inline
Rexp_inline=r;

options=optimset('Display','iter','TolX',1e-16,'display','off');
[y,~,~] =fzero(@ratio21,1,options);
if y >0
    BTd=y;
else
    y=log(P(nn:end));               % use log slope
    X=(1:length(y))/sampling_rate;
    BTd=exp(X/y);
end   
   
% given b calculate a from E1
e1=exp(1);
if (BTd==1)
    denom= (3-e1-1/e1);
else
    denom=(1-e1*exp(-BTd))/(BTd-1) - (e1-1)*(1-exp(-BTd))/BTd;
end
a=E1(end)/(Td*denom);

% Given b & a calculate c from E0
c = E0 - ((a./BTd) *(1 - exp(-BTd)));

% calculate b in s
b=BTd/Td;
pinf=c;
prd=a*exp(-b*td)+c;

% new fitting algorithm (31/07/08)
%options=optimset('TolX',1e-6);
%afind = @ (aa) sum((p(nn:end) - dias_int(p,t,aa,b,pinf,nn)).^2);
%aa=fminsearch(afind,0,options);
%pr=beat_int(p,t,aa,b,pinf,nn);
% alternative that works with earlier versions of Matlab
global p_inline
p_inline=p;
global t_inline
t_inline=t;
global b_inline
b_inline=b;
global pinf_inline
pinf_inline=pinf;
global nn_inline
nn_inline=nn;
global prd_inline
prd_inline=prd;

% options=optimset('TolX',1e-6);
options=optimset('TolX',1e-6,'display','off');
% following line can be used in versions that support explicit functions
%y=fminsearch(@dias_fit,[b],options);
[y,~,exitflag]=fminsearch(@dias_fit,b,options);
aa=y;
pr=beat_int(p,t,aa,b,pinf,nn);
if ~exitflag
    pr=min(p)*ones(1,length(p));
end

%find crossover point
Nd=length(pd);
np=[2:Nd Nd];
prend=pr(nn:end);
k=find((prend-prd).*(prend(np)-prd(np))<=0);
prr=pr;
prr(nn+k:end)=prd(k+1:end);

% return variables
Pr=prr;
A=aa ;
B=b;
Pinf=pinf;
Tn=t(nn);
Pn=p(nn);
end

%% ratio21.m - calculate the ratio of exponential moments R(bTd)
% KHP (27/01/10)
% added y==0 condition (06/02/10)

function Rdiff=ratio21(y)
%   input   y - value of bTd
%   output  R - ratio of E2/E1

global Rexp_inline

if (y==1)
    R=((exp(1)-1)-(1-exp(-1))*(exp(2)-1)/2)/(1-(1-exp(-1))*(exp(1)-1));
elseif (y==2)
    R=(1-(1-exp(-2))*(exp(2)-1)/4)/(-(exp(-1)-1)-(1-exp(-2))*(exp(1)-1)/2);
elseif (y==0)
    R=1/(3-exp(1));
else
    R=((exp(2-y)-1)/(2-y)-(1-exp(-y))*(exp(2)-1)/(2*y))/((exp(1-y)-1)/(1-y)-(1-exp(-y))*(exp(1)-1)/y);
end
Rdiff=R-Rexp_inline; 
end
 
%% kexpint.m  -  calculate the exponential integral Z = I_0^t Y*exp(at)dt
function Z=kexpint(Y,T,A)
%   inputs  Y(t) - row vector
%           T    - duration of Y
%           A    - exponential factor
%   Outputs Z    - row vector
N=length(Y)-1;
t=(0:N)*T/N;
dt=t(2)-t(1);
% calculate integral of y using trapezoidal rule
ye=Y.*exp(A*t);
Z=(cumtrapz(ye))*dt;
end

%% dias_int
% calculate the indefinite integral during whole beat, returning Pr_dias
% inputs    Ps   - row vector of pressure for whole beat
%           Ts   - row vector of times for whole beat
%           a    - fitting parameter
%           b    - reciprocal of the diastolic time constant tau
%           Pinf - asymptote of diastolic pressure
%           nn   - index of notch
% output    Prd  - Pr during diastole

function Prd = dias_int(Ps,Ts,a,b,Pinf,nn)

% modification to make fminsearch work in older versions of matlab
%function Pr = dias_int(a)
%global p_inline
%global t_inline
%global b_inline
%global Pinf_inline
%global nn_send_to_inline
%dt=t_inline(2)-t_inline(1);
%pse=cumtrapz(p_inline.*exp((a+b_inline)*ts_inline))*dt;
%prb=exp(-(a+b_inline)*t_inline).*(a*pse + p_inline(1) - (b_inline*Pinf_inline/(a+b_inline))) + b_inline*Pinf_inline/(a+b_inline);
%Pr=prb(nn:end);
%end

% calculate the integral using the trapezoidal rule
dt=Ts(2)-Ts(1);
pse=cumtrapz(Ps.*exp((a+b)*Ts))*dt;
prb = exp(-(a+b)*Ts).*(a*pse + Ps(1) - (b*Pinf/(a+b))) + b*Pinf/(a+b);
Prd=prb(nn:end);
end

%% beat_int
%  khp (30/07/08)
% calculate the indefinite integral during whole beat
% inputs    Ps   - row vector of pressure for whole beat
%           Ts   - row vector of times for whole beat
%           a    - fitting parameter
%           b    - reciprocal of the diastolic time constant tau
%           Pinf - asymptote of diastolic pressure
%           nn   - index of notch

function Pr = beat_int(Ps,Ts,a,b,Pinf,~)
% calculate the integral using the trapezoidal rule
dt=Ts(2)-Ts(1);
pse=cumtrapz(Ps.*exp((a+b)*Ts))*dt;
prb = exp(-(a+b)*Ts).*(a*pse + Ps(1) - (b*Pinf/(a+b))) + b*Pinf/(a+b);
Pr=prb;
end

%% dias_fit
function afind = dias_fit(y)

global p_inline
global t_inline
global b_inline
global pinf_inline
global nn_inline
global prd_inline

aa=y;
%afind = sum((p_inline(nn_inline:end) - dias_int(p_inline,t_inline,aa,b_inline,pinf_inline,nn_inline)).^2);
afind = sum((prd_inline - dias_int(p_inline,t_inline,aa,b_inline,pinf_inline,nn_inline)).^2);
end



%% Savitsky-Golay first derivative filter function
% dx=fsg721(x)
% 	7 point SavGol filter, 2nd order polynomial, 1st derivative
%	input 	x
%	output	dx
%	corrected for time shift

function dx=fsg721(x)

% 2nd order polynomial
C=[0.107143,0.071429,0.035714];

B=zeros(1,7);
for i=1:3 
    B(i)=C(i); 
end	
B(4)=0.0;
for i=5:7
    B(i)=-C(8-i);
end
A=[1,0];

s=size(x,2);
dx=filter(B,A,x);
dx=[dx(7),dx(7),dx(7),dx(7:s),dx(s),dx(s),dx(s)];

end
