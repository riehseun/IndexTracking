% MIE479 Capstone
% Seungmoon Rieh
% December, 2014
% 
% In the excel file, there is 501 stocks. Text array will consist of 501
% elements while ndata array will consist of 1506 elements. 
% test(1) to test(501) has strings representing each stock
% ndata(1) to ndata(501) has number representing sector that each stock 
% belongs.
% ndata(502) is undefined
% ndata(503) tp ndata(1003) has number representing market cap for each
% stock
% ndata(1004) is the total market cap for Index
% ndata(1005) to ndata(1505) has weight for each stock in Index
% ndata(1506) is the sum of weights 
filename = 'S&P500.xlsx';
[ndata, text] = xlsread(filename);

% S1 to S10 holds stock belonging to each industry sector
S1 = []; 
S2 = []; 
S3 = []; 
S4 = []; 
S5 = []; 
S6 = []; 
S7 = []; 
S8 = []; 
S9 = []; 
S10 = []; 

% s1 to s10 represent number of stocks to be included from each industry
% sector
s1 = 15;
s2 = 35;
s3 = 21;
s4 = 21;
s5 = 44;
s6 = 26;
s7 = 24;
s8 = 14;
s9 = 28;
s10 = 14;

% N = number of stocks in Index which is 501
N = size(text);
N = N(1); 

% Stratify stocks from each industry sector in order of market 
% capitalization 
i = 1; 
while i < N - 1
   %sector number 10: Communications 
   while ndata(i) == 10
       i = i + 1;
   end
   S1 = text(i-s1:i-1);
   %sector number 11: Comsumer discretionary
   while ndata(i) == 11
       i = i + 1;
   end
   S2 = text(i-s2:i-1);
   %sector number 12: Consumer staples
   while ndata(i) == 12
       i = i + 1;
   end
   S3 = text(i-s3:i-1);
   %sector number 13: Energy
   while ndata(i) == 13
       i = i + 1;
   end
   S4 = text(i-s4:i-1);
   %sector number 14: Financials
   while ndata(i) == 14
       i = i + 1;
   end
   S5 = text(i-s5:i-1);
   %sector number 15: Health care
   while ndata(i) == 15
       i = i + 1;
   end
   S6 = text(i-s6:i-1);
   %sector number 16: Industrials
   while ndata(i) == 16
       i = i + 1;
   end
   S7 = text(i-s7:i-1);
   %sector number 17: Materials
   while ndata(i) == 17
       i = i + 1;
   end
   S8 = text(i-s8:i-1);
   %sector number 18: Technology
   while ndata(i) == 18
       i = i + 1;
   end
   S9 = text(i-s9:i-1);
   %sector number 19: Utilities
   while ndata(i) == 19
       i = i + 1;
   end
   S10 = text(i-s10:i-1);
end

% Array to hold weights for stocks in S&P500
X = zeros(N,1);
for a = 2*N+3:3*N+2
    X(1+a-1005) = ndata(a); 
end

% Array to hold tickers for stock
stock = [S1;S2;S3;S4;S5;S6;S7;S8;S9;S10];

% generate comma separated list from cell array. This step is needed to 
% read prices from Yahoo Finance
n = s1+s2+s3+s4+s5+s6+s7+s8+s9+s10;
stocks = cell(1,n);
a = stock';
for i = 1:n
  stocks{i} = a{:,i};
end

% start date 
sD = '00';
sM = '01';
sY = '2004'; 

% end date
eD = '01';
eM = '01';
eY = '2014';

% using monthly prices
period = 'd';
dataFormat = 'cvs';

% number of periods
t = 2520;

% let p be the price matrix
pd = zeros(n,t+1);

for a = 1:n 
    urlname=['http://ichart.finance.yahoo.com/table.csv?'...
        '&s=' stocks{a} '&d=' eD '&e=' eM '&f=' eY '&g=' period '&a=' sD '&b=' sM '&c=' sY '&ignore=.' dataFormat];
    buff_reader=java.io.BufferedReader(java.io.InputStreamReader(openStream(java.net.URL(urlname))));
    % from the time period 1
    b = 1;
    % skip the first line
    line = buff_reader.readLine;
    while b<=t+1
        line = char(buff_reader.readLine);
        sepInd = strfind(line,',');
        % price of stock s at time t (most recent to past)
        try pd(a,b) = str2double(line(sepInd(6)+1:end)); catch end
        b = b+1;
    end
end

% save to file
save('dailyprice1.mat','pd','t','n','stocks','X');