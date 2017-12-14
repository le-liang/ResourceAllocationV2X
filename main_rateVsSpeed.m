% main function for D2D-based vehicular communications
% Compare sum/min ergodic capacity against varying vehicle velocity.

% By Le Liang, Georgia Tech, Jan. 26, 2017

tic
clear;
clc

channNum = 2e4;
rng(3); % control the random seed for randn, randi, rand

%% Parameters setup
infty = 2000; % used as infinity in the simulation
dB_Pd_max = 23; % max DUE transmit power in dBm
dB_Pc_max = 23; % max CUE transmit power in dBm

% large scale fading parameters
stdV2V = 3; % shadowing std deviation
stdV2I = 8;

% cell parameter setup
freq = 2; % carrier frequency 2 GHz
radius = 500; % cell radius in meters
bsHgt = 25; % BS height in meters
disBstoHwy = 35; % BS-highway distance in meters
bsAntGain = 8; % BS antenna gain 8 dBi
bsNoiseFigure = 5; % BS noise figure 5 dB

vehHgt = 1.5; % vehicle antenna height, in meters
vehAntGain = 3; % vehicle antenna gain 3 dBi
vehNoiseFigure = 9; % vehicle noise figure 9 dB

numLane = 6;
laneWidth = 4;
v = 60:10:140; % velocity
d_avg_ = 2.5.*v/3.6; % average inter-vehicle distance according to TR 36.885

% QoS parameters for CUE and DUE
r0 = 0.5; % min rate for CUE in bps/Hz
dB_gamma0 = 5; % SINR_min for DUE in dB
p0 = 0.001; % outage probability for DUE
dB_sig2 = -114; % noise power in dBm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% dB to linear scale conversion
sig2 = 10^(dB_sig2/10);
gamma0 = 10^(dB_gamma0/10);
Pd_max = 10^(dB_Pd_max/10);
Pc_max = 10^(dB_Pc_max/10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
numCUE = 20;
numDUE = 20;
sumRate_maxSum = zeros(length(d_avg_), 1);
minRate_maxSum = zeros(length(d_avg_), 1);
sumRate_maxMin = zeros(length(d_avg_), 1);
minRate_maxMin = zeros(length(d_avg_), 1);
%%

parfor ind = 1 : length(sumRate_maxSum)
    d_avg = d_avg_(ind);
    
    cntChann = 0; % channel realization counter
    while cntChann < channNum
        %% Generate traffic on the highway
        d0 = sqrt(radius^2-disBstoHwy^2);
        [genFlag,vehPos,indCUE,indDUE,indDUE2] = genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE, numDUE);
        if genFlag == 1
            continue; % generated vehicles are not enough to perform simulation, jump to the next iteration.
        end
        
        %% random large-scale fading generation
        alpha_mB_ = zeros(1, numCUE);
        alpha_k_ = zeros(1, numDUE);
        alpha_kB_ = zeros(1, numDUE);
        alpha_mk_ = zeros(numCUE, numDUE);
        for m = 1 : numCUE
            dist_mB = sqrt(vehPos(indCUE(m),1)^2 + vehPos(indCUE(m),2)^2);
            dB_alpha_mB = genPL('V2I', stdV2I, dist_mB, vehHgt, bsHgt, freq) + vehAntGain+bsAntGain-bsNoiseFigure;
            alpha_mB_(m) = 10^(dB_alpha_mB/10);
            
            for k = 1 : numDUE
                dist_mk = sqrt((vehPos(indCUE(m),1)-vehPos(indDUE(k),1))^2 +  (vehPos(indCUE(m),2)-vehPos(indDUE(k),2))^2);
                dB_alpha_mk = genPL('V2V', stdV2V, dist_mk, vehHgt, vehHgt, freq) + 2*vehAntGain-vehNoiseFigure;
                alpha_mk_(m,k) = 10^(dB_alpha_mk/10);
            end
        end
        for k = 1 : numDUE
            dist_k = sqrt((vehPos(indDUE(k),1)-vehPos(indDUE2(k),1))^2 +  (vehPos(indDUE(k),2)-vehPos(indDUE2(k),2))^2);
            dB_alpha_k = genPL('V2V', stdV2V, dist_k, vehHgt, vehHgt, freq) + 2*vehAntGain-vehNoiseFigure;
            alpha_k_(k) = 10^(dB_alpha_k/10);
            
            dist_k = sqrt(vehPos(indDUE(k),1)^2 + vehPos(indDUE(k),2)^2);
            dB_alpha_kB = genPL('V2I', stdV2I, dist_k, vehHgt, bsHgt, freq)+ vehAntGain+bsAntGain-bsNoiseFigure;
            alpha_kB_(k) = 10^(dB_alpha_kB/10);
        end
        
        %% resource allocation design - single pair
        C_mk = zeros(numCUE, numCUE); % create virtual DUEs if numDUE < numCUE
        for m = 1 : numCUE
            alpha_mB = alpha_mB_(m);
            for k = 1 : numCUE
                
                if k > numDUE % create virtual DUEs if numDUE < numCUE
                    a = Pc_max*alpha_mB/sig2;% in absence of DUEs, the CUEs transmit at max power
                    C_mk(m,k) = computeCapacity(a,0); % no interference from DUEs
                    continue;
                end
                
                alpha_k = alpha_k_(k);
                alpha_kB = alpha_kB_(k);
                alpha_mk = alpha_mk_(m,k);
                
                Pc_dmax = alpha_k*Pd_max/(gamma0*alpha_mk)*(exp(-gamma0*sig2/(Pd_max*alpha_k))/(1-p0)-1);
                if Pc_dmax <= Pc_max
                    Pd_opt = Pd_max;
                    Pc_opt = Pc_dmax;
                else
                    %% Bisectin search to find Pd_cmax
                    epsi = 1e-5;
                    Pd_left = -gamma0*sig2/(alpha_k*log(1-p0)); % P_{k,min}^d
                    Pd_right = Pd_max;
                    tmpVeh = 0;
                    while Pd_right - Pd_left > epsi
                        tmpVeh = (Pd_left + Pd_right)/2;
                        if alpha_k*tmpVeh/(gamma0*alpha_mk)*(exp(-gamma0*sig2/(tmpVeh*alpha_k))/(1-p0)-1) > Pc_max
                            Pd_right = tmpVeh;
                        else
                            Pd_left = tmpVeh;
                        end
                    end
                    Pd_cmax = tmpVeh;
                    
                    Pd_opt = Pd_cmax;
                    Pc_opt = Pc_max;
                end
                %% C_km stores the optimal throughput for the "m"th CUE when sharing
                % spectrum with the "k"th DUE
                a = Pc_opt*alpha_mB/sig2;
                b = Pd_opt*alpha_kB/sig2;
                C_mk(m,k) = computeCapacity(a,b);
                if C_mk(m,k) < r0 % min rate for the V2I link
                    C_mk(m,k) = -infty;
                end
            end
        end
        
        %% Reuse pair matching
        [assignmentSum, cost] = munkres(-C_mk);
        [sumVal_sum, minVal_sum] = sumAndMin(C_mk, assignmentSum);
        
        [assignmentMin, dummyMin ] = maxMin( C_mk );
        [sumVal_min, minVal_min] = sumAndMin(C_mk, assignmentMin);
        
        if minVal_sum < 0 ||  minVal_min < 0 % infeasible problem
            continue;
        end
        
        sumRate_maxSum(ind) = sumRate_maxSum(ind) + sumVal_sum;
        minRate_maxSum(ind) = minRate_maxSum(ind) + minVal_sum;
        sumRate_maxMin(ind) = sumRate_maxMin(ind) + sumVal_min;
        minRate_maxMin(ind) = minRate_maxMin(ind) + minVal_min;
    
        cntChann = cntChann + 1;
    end
    ind
    
end

sumRate_maxSum = sumRate_maxSum/channNum;
minRate_maxSum = minRate_maxSum/channNum;
sumRate_maxMin = sumRate_maxMin/channNum;
minRate_maxMin = minRate_maxMin/channNum;
%%



%% Second run with different max powers
dB_Pd_max = 17; % max DUE transmit power in dBm
dB_Pc_max = 17; % max CUE transmit power in dBm
Pd_max = 10^(dB_Pd_max/10);
Pc_max = 10^(dB_Pc_max/10);

sumRate_maxSum2 = zeros(length(d_avg_), 1);
minRate_maxSum2 = zeros(length(d_avg_), 1);
sumRate_maxMin2 = zeros(length(d_avg_), 1);
minRate_maxMin2 = zeros(length(d_avg_), 1);
%%
parfor ind = 1 : length(sumRate_maxSum)
    d_avg = d_avg_(ind);
    
    cntChann = 0; % channel realization counter
    while cntChann < channNum
        %% Generate traffic on the highway
        d0 = sqrt(radius^2-disBstoHwy^2);
        [genFlag,vehPos,indCUE,indDUE,indDUE2] = genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE,numDUE);
        if genFlag == 1
            continue; % generated vehicles are not enough to perform simulation, jump to the next iteration.
        end
        
        %% random large-scale fading generation
        alpha_mB_ = zeros(1, numCUE);
        alpha_k_ = zeros(1, numDUE);
        alpha_kB_ = zeros(1, numDUE);
        alpha_mk_ = zeros(numCUE, numDUE);
        for m = 1 : numCUE
            dist_mB = sqrt(vehPos(indCUE(m),1)^2 + vehPos(indCUE(m),2)^2);
            dB_alpha_mB = genPL('V2I', stdV2I, dist_mB, vehHgt, bsHgt, freq) + vehAntGain+bsAntGain-bsNoiseFigure;
            alpha_mB_(m) = 10^(dB_alpha_mB/10);
            
            for k = 1 : numDUE
                dist_mk = sqrt((vehPos(indCUE(m),1)-vehPos(indDUE(k),1))^2 +  (vehPos(indCUE(m),2)-vehPos(indDUE(k),2))^2);
                dB_alpha_mk = genPL('V2V', stdV2V, dist_mk, vehHgt, vehHgt, freq) + 2*vehAntGain-vehNoiseFigure;
                alpha_mk_(m,k) = 10^(dB_alpha_mk/10);
            end
        end
        for k = 1 : numDUE
            dist_k = sqrt((vehPos(indDUE(k),1)-vehPos(indDUE2(k),1))^2 +  (vehPos(indDUE(k),2)-vehPos(indDUE2(k),2))^2);
            dB_alpha_k = genPL('V2V', stdV2V, dist_k, vehHgt, vehHgt, freq) + 2*vehAntGain-vehNoiseFigure;
            alpha_k_(k) = 10^(dB_alpha_k/10);
            
            dist_k = sqrt(vehPos(indDUE(k),1)^2 + vehPos(indDUE(k),2)^2);
            dB_alpha_kB = genPL('V2I', stdV2I, dist_k, vehHgt, bsHgt, freq)+ vehAntGain+bsAntGain-bsNoiseFigure;
            alpha_kB_(k) = 10^(dB_alpha_kB/10);
        end
        
        %% resource allocation design - single pair
        C_mk = zeros(numCUE, numCUE); % create virtual DUEs if numDUE < numCUE
        for m = 1 : numCUE
            alpha_mB = alpha_mB_(m);
            for k = 1 : numCUE
                
                if k > numDUE % create virtual DUEs if numDUE < numCUE
                    a = Pc_max*alpha_mB/sig2;% in absence of DUEs, the CUEs transmit at max power
                    C_mk(m,k) = computeCapacity(a,0); % no interference from DUEs
                    continue;
                end
                
                alpha_k = alpha_k_(k);
                alpha_kB = alpha_kB_(k);
                alpha_mk = alpha_mk_(m,k);
                
                Pc_dmax = alpha_k*Pd_max/(gamma0*alpha_mk)*(exp(-gamma0*sig2/(Pd_max*alpha_k))/(1-p0)-1);
                if Pc_dmax <= Pc_max
                    Pd_opt = Pd_max;
                    Pc_opt = Pc_dmax;
                else
                    %% Bisectin search to find Pd_cmax
                    epsi = 1e-5;
                    Pd_left = -gamma0*sig2/(alpha_k*log(1-p0)); % P_{k,min}^d
                    Pd_right = Pd_max;
                    tmpVeh = 0;
                    while Pd_right - Pd_left > epsi
                        tmpVeh = (Pd_left + Pd_right)/2;
                        if alpha_k*tmpVeh/(gamma0*alpha_mk)*(exp(-gamma0*sig2/(tmpVeh*alpha_k))/(1-p0)-1) > Pc_max
                            Pd_right = tmpVeh;
                        else
                            Pd_left = tmpVeh;
                        end
                    end
                    Pd_cmax = tmpVeh;
                    
                    Pd_opt = Pd_cmax;
                    Pc_opt = Pc_max;
                end
                %% C_km stores the optimal throughput for the "m"th CUE when sharing
                % spectrum with the "k"th DUE
                a = Pc_opt*alpha_mB/sig2;
                b = Pd_opt*alpha_kB/sig2;
                C_mk(m,k) = computeCapacity(a,b);
                if C_mk(m,k) < r0 % min rate for the V2I link
                    C_mk(m,k) = -infty;
                end
            end
        end
        
        %% Reuse pair matching
        [assignmentSum, cost] = munkres(-C_mk);
        [sumVal_sum, minVal_sum] = sumAndMin(C_mk, assignmentSum);
        
        [assignmentMin, dummyMin ] = maxMin( C_mk );
        [sumVal_min, minVal_min] = sumAndMin(C_mk, assignmentMin);
        
        if minVal_sum < 0 ||  minVal_min < 0 % infeasible problem
            continue;
        end
        
        sumRate_maxSum2(ind) = sumRate_maxSum2(ind) + sumVal_sum;
        minRate_maxSum2(ind) = minRate_maxSum2(ind) + minVal_sum;
        sumRate_maxMin2(ind) = sumRate_maxMin2(ind) + sumVal_min;
        minRate_maxMin2(ind) = minRate_maxMin2(ind) + minVal_min;
    
        cntChann = cntChann + 1;
    end
    ind
end

sumRate_maxSum2 = sumRate_maxSum2/channNum;
minRate_maxSum2 = minRate_maxSum2/channNum;
sumRate_maxMin2 = sumRate_maxMin2/channNum;
minRate_maxMin2 = minRate_maxMin2/channNum;


%%
LineWidth = 1.5;
MarkerSize = 6;
figure
plot(v, sumRate_maxSum, 'k-s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, sumRate_maxMin, 'b-o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, sumRate_maxSum2, 'k--s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, sumRate_maxMin2, 'b--o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
grid on
legend('P_{max}^c = 23 dBm, Algorithm 1', 'P_{max}^c = 23 dBm, Algorithm 2',...
'P_{max}^c = 17 dBm, Algorithm 1', 'P_{max}^c = 17 dBm, Algorithm 2')
xlabel('$v$ (km/h)', 'interpreter','latex')
ylabel('$\sum\limits_m C_m$ (bps/Hz)', 'interpreter','latex')
% saveas(gcf, sprintf('sumRateVsSpeed')); % save current figure to file

figure
plot(v, minRate_maxSum, 'k-s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, minRate_maxMin, 'b-o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, minRate_maxSum2, 'k--s', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
hold on
plot(v, minRate_maxMin2, 'b--o', 'LineWidth', LineWidth, 'MarkerSize', MarkerSize)
grid on
legend('P_{max}^c = 23 dBm, Algorithm 1', 'P_{max}^c = 23 dBm, Algorithm 2',...
'P_{max}^c = 17 dBm, Algorithm 1', 'P_{max}^c = 17 dBm, Algorithm 2')
xlabel('$v$ (km/h)', 'interpreter','latex')
ylabel('$\min C_m$ (bps/Hz)', 'interpreter','latex')
% saveas(gcf, 'minRateVsSpeed'); % save current figure to file

% save all data
save('main_rateSpeed_Jan29')

toc