% function [H] = channel(param)
% 
% R = param.no_antenna;
% L = param.no_carrier;

R = 100;     % number of antennas 
L = 1024;    % number of subcariers

s = qd_simulation_parameters; 
s.center_frequency = [3.5e9];    % central frequency
i = qd_layout( s ); 
i.tx_position = [0 0 25]';       % position of transmitter, the height of transmitter is 25m
i.no_rx = 1;                     % number of receiver   
i.rx_position = [200 40 1.5]';   % position of receiver
indoor_rx = i.set_scenario('3GPP_38.901_UMa_NLOS',[],[],0.8);   % use "3GPP_38.901_UMa_NLOS" setup
a_35000_MHz = qd_arrayant( '3gpp-3d', R, 1, s.center_frequency(1), 1 );
s.show_progress_bars = 0;    % Disable progress bar
i.tx_array(1,1) = a_35000_MHz;
i.rx_array = qd_arrayant('omni');

sample_distance = 5;
x_min = -50;
x_max = 550;
y_min = -300;
y_max = 300;
rx_height = 1.5;             % height of receiver

s.show_progress_bars = 0;    % Disable progress bar
c = i.get_channels;          % generate the time-domain channel 

% transform the channel from time domain to frequency domain 
B  = 100e6;  % bandwith 10MHz

H_freq = c(1,1).fr(B,L);  % transform channel from time domain to frequency domain
[idx1,idx2,idx3] = size(H_freq);
H_af = zeros(idx2,idx3);
H_af(:,:) = H_freq(1,:,:);   % Get the antenna and frequency domain channel matrix with dimension R-by-L

% define DFT matrix
% antenna domain to angular domain DFT 
FB = zeros(R,R);
for i = 1:R
    for j = 1:R
        FB(i,j) = 1/sqrt(R)*exp(1)^(2*pi*(i-1)*(j-1)*sqrt(-1)/R);
    end
end

% frequency domain to delay domain DFT
FC = zeros(L,L);
for i = 1:L
    for j = 1:L
        FC(i,j) = 1/sqrt(L)*exp(1)^(2*pi*(i-1)*(j-1)*sqrt(-1)/L);
    end 
end


% Get the antenna and delay domain channel matrix
H_ad = H_af*FC;


% generate the channel sparsity 
% sc=1 sp=8 si=1

% generate partial sparsity sp = 8
H_ad(1:R/4*3,2:9)     = 0;
H_ad(1:R/2,10:17)     = 0;
H_ad(R/4*3+1:R,10:17) = 0;
H_ad(1:R/4,18:25)     = 0;
H_ad(R/2+1:R,18:25)   = 0;
H_ad(R/4+1:R,26:33)   = 0;

% generate individual sparsity si = 1
H_ad(:,34:33+R) = H_ad(:,34:33+R).*eye(R);

% set the remaining elements to 0
H_ad(:,34+R:L) = 0;





