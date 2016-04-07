% MIE479 Capstone
% Seungmoon Rieh
% December, 2014

% load simulated wAVaR of Index
load wAVaR.mat;
u_estimate = mean(wAVaR);
s = std(wAVaR);
% 90% z score from normal dist = 1.645
xp = u_estimate + 1.645*s/sqrt(100000);
xn = u_estimate - 1.645*s/sqrt(100000);
u_hat = xp - u_estimate;
% save to file
save('ci.mat','u_estimate','u_hat');

