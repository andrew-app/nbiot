close all

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


rmit = [-37.808176661814905, 144.96240622449]; %location of building 80


eff_nf = -85;

SNR_min_c1 = 6.28; 
%from BLER vs SNR graph 
%where BLER is 5% for 1 rr

SNR_min_c2 = -10; 
%from BLER vs SNR graph
%where BLER is 5% for 32 rr

rx_sensitivity_c1 = eff_nf + SNR_min_c1; 

%Receiver sensitivity for 1 rr

rx_sensitivity_c2 = eff_nf + SNR_min_c2; 

%Receiver sensitivity for 32 rr

SNR_min = -10; %from BLER vs SNR graph for both repetitons SNR where BLER is 5%





rx = rxsite('Name','RMIT Building 80', ...
       'Latitude', rmit(1), ...
       'Longitude', rmit(2), ...
       'ReceiverSensitivity', rx_sensitivity_c2); %Usage of 32 rr since better results from simulation

ss = sigstrength(rx,txs);

margin = abs(rx_sensitivity_c2 - ss);


figure(1)

pattern(array,fc)

sites = 1:length(ss);

figure(2)

plot(sites,ss);
xlabel('Site Number')
ylabel('RSS(dBm)')




distgc = zeros(1,length(ss));
dist = zeros(1,length(ss));
max_rss = max(ss);
min_rss = min(ss);


for j = 1:length(ss)
    
    if ss(j) == max_rss
        max_rss = [td(j,1),td(j,2)];
        best_site = j;
    elseif ss(j) == min_rss
        worst_site = j;    
        min_rss = [td(j,1),td(j,2)];
    end
    
    
end

for j = 1:length(ss)
    
    distgc(j) = distance(td(j,1),td(j,2),rmit(1),rmit(2));
    
    
    dist(j) = deg2km(distgc(j));
    
    
    
    
end



figure(3)
geobubble(td(:,1),td(:,2))



rx_sites = [max_rss; min_rss; rmit];


figure(4)
geobubble(rx_sites(:,1),rx_sites(:,2)) %map showing location of lowest and highest rssi tx relative to receiver






dist = transpose(dist);

X = [ones(length(dist),1) dist];


m = X\ss;

y = X*m;
figure(5)
scatter(dist,ss)
xlabel('distance(km)')
ylabel('RSSI (dBm)')
hold on
plot(dist,y)
hold off




%coverage(txs, 'longley-rice', 'SignalStrengths',-110:10:0) %RSS Map

%sinr(txs, 'Values',[-10,10]) %SNR MAP
%coverage(txs, 'longley-rice','Resolution', 350, 'MaxRange', 6000, 'SignalStrengths', rx_sensitivity_c1, 'Colors', 'blue', 'ReceiverGain', 0, 'ReceiverAntennaHeight', 2) %1rr
%coverage(txs, 'longley-rice','Resolution', 350, 'MaxRange',%6000,'SignalStrengths', rx_sensitivity_c2, 'ReceiverGain', 0,'ReceiverAntennaHeight', 2) 
%32rr