%modelfun = @(b,x)(b(1)+exp(-b(2)*x));
modelfun = @(b,t)(exp(-b(1)*t)+b(2));
rng('default') % for reproducibility
b = [1;1];
t = exprnd(2,100,1);
y = modelfun(b,t) + normrnd(0,0.1,100,1);
y= y+20;
opts.RobustWgtFun = 'bisquare';
opts1.RobustWgtFun = 'huber';
opts2.RobustWgtFun = 'welsch'
plot(t,y,'.')
beta0 = [2;20];
beta = nlinfit(t,y,modelfun,beta0)          % none
beta1 = nlinfit(t,y,modelfun,beta0,opts)    % bisquare
beta2 = nlinfit(t,y,modelfun,beta0,opts1)   % huber
beta3 = nlinfit(t,y,modelfun,beta0,opts2)   % welsch

% Huber looks best on this test