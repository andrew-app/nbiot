td = csvread('TowersLocations.csv');

load('Directional12dBi.mat')


loc = geobubble(td(:,1),td(:,2));

tx = txsite(  ...
    'Latitude', td(:,1), ...
    'Longitude', td(:,2), ...
    'AntennaHeight', 30, ...
    'TransmitterPower', 40, ...
    'TransmitterFrequency', 900.0e6 ...
);





coverage(tx, 'longley-rice', 'SignalStrengths', -80);
