% MIE479 Capstone
% Seungmoon Rieh
% December, 2014

% load X, weights of assets in Index
load dailyprice1.mat;
j = 10;

% calculate weights for observation at tail
w = zeros(j,1);
ARA = 1;
for a = 1:j
    w(a,1) = ARA*exp(ARA*(1-((a-1)/j)))/(exp(ARA)-1) - ARA*exp(ARA*(1-a/j))/(exp(ARA)-1);
end

wAVaR = zeros(10,1);
for c = 1:100000
    k = size(X);
    I = k(1); % number of assets
    J = 10*j; % total number of observation for 90% interval
    r = trnd(3,I,J)/10; % generate random daily return for each stock 
    rp = (X'*r); % return for Index
    rp = sort(rp,'ascend'); % sort return such that most negative is at top
    a = 0.9;
    level = round((1-a)*J); % confidence level 
    rp = rp(:,1:level); % pick return at tail
    VaR = rp(level); % value at risk of return for Index
    w = sort(w,'descend'); % sort weight such that highest is at top
    wAVaR(c) = rp*w; % weighted average value of return for Index
end

% save to file
save('wAVaR.mat','wAVaR');