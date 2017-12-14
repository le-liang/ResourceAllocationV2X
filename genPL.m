function [ combinedPL ] = genPL(linkType, stdShadow, dist, hgtTX, hgtRX, freq)
% genPL: compute large scale fading (pathloss + shadowing) parameters
% Input: mode = 'v2v' or 'v2i', specify link types of vehcular commun.
%        stdShadow, log-normal shadowing std deviation in dB
%        dist, distance between TX and RX in meters
%        heightTX, transmitter antenna height in meters
%        heightRX, receiver antenna height in meters
%        freq, carrier frequency in GHz
% Output: combinedPL, combined large scale fading value (pathloss +
% shadow) in dB
% By Le Liang, Georgia Tech, July 10, 2016

if strcmp(upper(linkType), 'V2V')
    d_bp = 4*(hgtTX-1)*(hgtRX-1)*freq*10^9/(3*10^8);
    A = 22.7; B = 41.0; C = 20;
    if dist <= 3
        PL = A*log10(3) + B + C*log10(freq/5);
    elseif dist <= d_bp
        PL = A*log10(dist) + B + C*log10(freq/5);
    else
        PL = 40*log10(dist)+9.45-17.3*log10((hgtTX-1)*(hgtRX-1))+2.7*log10(freq/5);
    end
else
    PL = 128.1 + 37.6*log10(sqrt((hgtTX-hgtRX)^2+dist^2)/1000);
end

combinedPL = - (randn(1)*stdShadow + PL);

end