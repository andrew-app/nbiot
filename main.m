td = csvread('TowersLocations.csv');

load('Directional12dBi.mat')





txs = txsite(  ...
    'Latitude', td(:,1), ...
    'Longitude', td(:,2), ...
    'AntennaHeight',30,...
    'TransmitterPower', 40, ...
    'TransmitterFrequency', fc);



array = phased.UCA('Element',antenna,'Radius',2,'NumElements',3);


for i = 1:length(txs)
   txs(i).Antenna = array; 
end

bw = 180000;


rmit = [-37.808176661814905, 144.96240622449];


eff_nf = -85;

SNR_min = -22; %from BLER vs SNR graph for 32 repetitons first SNR where BLER is 100%

rx_sensitivity = eff_nf + SNR_min;



rx = rxsite('Name','RMIT Building 80', ...
       'Latitude', rmit(1), ...
       'Longitude', rmit(2), ...
       'ReceiverSensitivity', rx_sensitivity);

ss = sigstrength(rx,txs);

margin = abs(rx_sensitivity - ss);


figure(1)

pattern(array,fc)

sites = 1:length(ss);

figure(2)

plot(sites,ss);
xlabel('Site Number')
ylabel('RSS(dBm)')




distgc = zeros(1,length(ss));
distrh = zeros(1,length(ss));
dist = zeros(1,length(ss));
max_rss = max(ss);
min_rss = min(ss);
k = 0;
for j = 1:length(ss)
    
    if ss(j) == max_rss
        max_rss = [td(j,1),td(j,2)];
    elseif ss(j) == min_rss
        min_rss = [td(j,1),td(j,2)];
    end
    
    
end
figure(2)
geobubble(td(:,1),td(:,2))



rx_sites = [max_rss; min_rss; rmit];

figure(3)
geobubble(rx_sites(:,1),rx_sites(:,2))
%coverage(txs,rx, 'longley-rice','SignalStrengths',-110:10:0)
%sinr(txs,)
