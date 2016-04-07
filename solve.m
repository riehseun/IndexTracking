% MIE479 Capstone
% Seungmoon Rieh
% December, 2014

% load historical price of selected stokcs, SPDR, and S%P500
load dailyprice1.mat;
load dailyprice2.mat;
% load wAVaR estimate and range
load ci.mat;

% Grap one year data of daily returns 
% Change this to test for different time periods
p = pd(:,1:253);
pp = ppd(:,1:253);

% Grap time period (number of trading days)
t = 252;
tt = 252;

% Grap wAVaR estimate and range
wAVaR = u_estimate;
l = u_hat;

% Calculating geometrical average daily return for stocks. 
% For stocks without historical price information, set u = 0
u = zeros(n,1);
for a = 1:n
   if p(a,t) == 0
       u(a) = 0;
   else
       u(a) = (p(a,1) - p(a,t)) / p(a,t);  
   end
end
u = (u+1).^(1/t)-1;

% Some stocks do not have historical price information. 
% Let index matrix contain stocks without historical
% price information bewteen chosen time period
index = [];
i = 1;
for a = 1:n
    if u(a) == 0
        index(i) = a;
        i = i + 1;
    end
end

% Calculating geometrical average daily return for stocks after removing
% stocks without historical price information. 
u = removerows(u,'ind',index);

% Calculating historical daily return for 252 days
e = size(index);
e = e(2);
r = zeros(n-e,t);
for a = 1:n 
   for b = 1:t
       r(a,b) = ((p(a,b) - p(a,b+1))/p(a,b));
   end
end

% remove stocks without historical price information
r = removerows(r,'ind',index);

% Stock tickers 
tickers = stocks';
tickers = removerows(tickers,'ind',index);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i and j will completely determine the number of variables and constraints
% in optimization framework

% Number of stocks in portfolio
i = size(u); 
i = i(1);

% Number of observations at tail.
% This number is most important input parameter. Large number correspond 
% to more observation at tail, which leads to better accuracy. However, 
% increase in j will drastically increase the number of constraints in 
% linprog, which may lead MATLAB to run out of memory
j = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimate average daily return of each stock
v_gen = trnd(3,i,100000)/10;
v = mean(v_gen,2);     

% Calculate weights for portfolio return observation at tail
w = zeros(j,1);
ARA = 1;
for a = 1:j
    w(a,1) = ARA*exp(ARA*(1-((a-1)/j)))/(exp(ARA)-1) - ARA*exp(ARA*(1-a/j))/(exp(ARA)-1);
end

% Estimate average daily return of each stock for j observations at tail
rij = zeros(i,j);
rij_hat = zeros(i*j,1);
for a = 1:j
    r_gen = trnd(3,i,100000)/10;
    r_estimate = mean(r_gen,2);     
    s = std(r_gen,0,2);
    % 90% z score from normal dist = 1.645
    xp = r_estimate + 1.645*s/sqrt(100000);
    xn = r_estimate - 1.645*s/sqrt(100000);
    rij(:,a) = [r_estimate];
    rij_hat(i*(j-1)+1:i*j) = [xp - r_estimate]/10;
end

% Input parameters for optimization framework
K = 0.05;
gamma = 30*ones(j,1);
c = 1;

% Objective function
f = [zeros(i,1);zeros(j,1);zeros(i*j,1);zeros(j,1);zeros(i,1);zeros(j,1);0;1]';

% Inequality constraints
A1 = [-v;zeros(j,1);zeros(i*j,1);zeros(j,1);zeros(i,1);zeros(j,1);0;-1]';
A2 = [zeros(i,1);-ones(j,1);zeros(i*j,1);zeros(j,1);zeros(i,1);w;K+1;0]';
A3 = [zeros(i,1);ones(j,1);zeros(i*j,1);zeros(j,1);zeros(i,1);-w;K-1;0]';

A4 = zeros(j,i+j+i*j+j+i+j+1+1);
W = diag(w);
for a = 1:j
    D = diag(w(a)*ones(i,1));
    A4(a,:) = [D*rij(:,a);-ones(j,1);-ones(i*j,1);-gamma;zeros(i,1);-W(:,a)+zeros(j,1);0;0];
end

A5 = zeros(j-1,i+j+i*j+j+i+j+1+1);
for a = 1:j-1
    A5(a,:) = [rij(:,a);zeros(j,1);-ones(i*j,1);-gamma;zeros(i,1);0;ones(j-1,1)*w(a);0;0];
end

A6 = zeros(j-1,i+j+i*j+j+i+j+1+1);
cof = 0.8;
mat1 = diag(ones(j,1));
mat2 = diag(cof*ones(j-1,1));
mat1(:,j) = [];
mat2 = [zeros(j-1,1)';mat2];
for a = 1:j-1
    A6(a,:) = [zeros(i,1);zeros(j,1);zeros(i*j,1);zeros(j,1);zeros(i,1);mat1(:,a)-mat2(:,a);0;0];
end

A7 = zeros(i*j,i+j+i*j+j+i+j+1+1);
for a = 1:i*j
    A7(a,:) = [zeros(i,1);zeros(j,1);-ones(i*j,1);-ones(j,1);rij_hat(a)*ones(i,1);zeros(j,1);0;0];
end

A8 = zeros(i,i+j+i*j+j+i+j+1+1);
D = diag(ones(i,1));
for a = 1:i
    A8(a,:) = [D(:,a);zeros(j,1);zeros(i*j,1);zeros(j,1);-D(:,a);zeros(j,1);0;0];
end

A9 = zeros(i,i+j+i*j+j+i+j+1+1);
for a = 1:i
    A9(a,:) = [-D(:,a);zeros(j,1);zeros(i*j,1);zeros(j,1);-D(:,a);zeros(j,1);0;0];
end

A = [A1;A2;A3;A4;A5;A6;A7;A8;A9];
b = [0;(K+1)*wAVaR;(K-1)*wAVaR;zeros(j,1);zeros(j-1,1);zeros(j-1,1);zeros(i*j,1);zeros(i,1);zeros(i,1)];

% Equality constraints
Aeq = [ones(i,1);zeros(j,1);zeros(i*j,1);zeros(j,1);zeros(i,1);zeros(j,1);0;0]';
beq = 1;

% lower and upper bounds
lb = [0.0002*ones(i,1);zeros(j,1);zeros(i*j,1);zeros(j,1);-Inf(i,1);-c;-Inf(j-1,1);-l;-Inf];
ub = [0.02*ones(i,1);Inf(j,1);Inf(i*j,1);Inf(j,1);Inf(i,1);c;Inf(j-1,1);l;Inf];

% Solve
options = optimset('LargeScale','off','Simplex','on');
x0 = 0.0002*ones(i+j+i*j+j+i+j+1+1,1);
[x,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub,x0,options);
res = x(1:i,1);

% Compute return of portfolio
rp = u'*res;

% Compute tail observations of portfolio 
pr = res'*r;
pr = sort(pr,'ascend'); % most negative is at top
len = size(pr); 
len = len(2);
lev = round(0.9*t);
prr = pr(:,1:t-lev); % pick tail observation

% calculate weights 
j= t-lev;
wp = zeros(j,1);
ARA = 1;
for a = 1:j
    wp(a,1) = ARA*exp(ARA*(1-((a-1)/j)))/(exp(ARA)-1) - ARA*exp(ARA*(1-a/j))/(exp(ARA)-1);
end

% Compute wAVaR of portfolio
wavarp = prr*wp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute SPDR and S&P500

% daily returns for SPDR and S&P500 for 252 days
rr = zeros(nn,tt);
for a = 1:nn
   for b = 1:tt
       rr(a,b) = ((pp(a,b) - pp(a,b+1))/pp(a,b));
   end
end 

% geometrical average daily return
uu = zeros(nn);
for a = 1:nn   
    uu(a) = (pp(a,1) - pp(a,tt)) / pp(a,tt);  
end
uu = prod(uu+1,2).^(1/tt)-1;

% confidence level
a = 0.9;
level = round((1-a)*tt);
rr = sort(rr,2,'ascend'); % sort such that most negative is at top
rr = rr(:,1:level); % pick observations at tail

% compute weights (highest to lowest)
ww = zeros(level,1);
ARA = 1;
for a = 1:level
    ww(a,1) = ARA*exp(ARA*(1-((a-1)/level)))/(exp(ARA)-1) - ARA*exp(ARA*(1-a/level))/(exp(ARA)-1);
end
ww = diag(ww);

z = rr*ww; %first col = SPDR, second col = SP500
wavar = sum(z,2);
rs = uu(1); %return of SPDR
ri = uu(2); %return of S&P500
wavars = wavar(1); %wAVaR of SPDR
wavari = wavar(2); %wAVaR of S&P500
TRs = abs((wavar(1)-wavar(2))/wavar(2)); %tracking error of SPDR
TRp = abs((wavarp-wavar(2))/wavar(2)); %tracking error of portfolio

[rs rp ri wavars wavarp wavari TRs TRp]