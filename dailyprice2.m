% MIE479 Capstone
% Seungmoon Rieh
% December, 2014

% SPDR ETF, S&P500 Index 
etf = {'SPY','^GSPC'};

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

% number of stocks
d = size(etf);
nn = d(2);

% number of periods 16*12 = 192 / 16*252 = 4032
tt = 2520;

% let p be the price matrix
ppd = zeros(nn,tt+1);

for a = 1:nn 
    urlname=['http://ichart.finance.yahoo.com/table.csv?'...
        '&s=' etf{a} '&d=' eD '&e=' eM '&f=' eY '&g=' period '&a=' sD '&b=' sM '&c=' sY '&ignore=.' dataFormat];
    buff_reader=java.io.BufferedReader(java.io.InputStreamReader(openStream(java.net.URL(urlname))));
    % from the time period 1
    b = 1;
    % skip the first line
    line = buff_reader.readLine;
    while b<=tt+1
        line = char(buff_reader.readLine);
        sepInd = strfind(line,',');
        % price of stock s at time t (most recent to past)
        try ppd(a,b) = str2double(line(sepInd(6)+1:end)); catch end
        b = b+1;
    end
end

% save to file
save('dailyprice2.mat','ppd','tt','nn','etf');