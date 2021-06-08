td = csvread('TowersLocations.csv');

load('Directional12dBi.mat')





tx = txsite(  ...
    'Latitude', td(:,1), ...
    'Longitude', td(:,2), ...
    'TransmitterPower', 40, ...
    'TransmitterFrequency', fc);



array = phased.UCA('Element',antenna,'Radius',2,'NumElements',3);


for i = 1:length(tx)
   tx(i).Antenna = array; 
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

ss = sigstrength(rx,tx);

margin = abs(rx_sensitivity - ss);


figure(1)

pattern(array,fc)


sinr(tx)

distgc = distance(-15,0,60,150);

