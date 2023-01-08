% Niall Donohoe
% 17340571
% Digital Communications - Assignment 2
% Simulation of a Hamming-coded BPSK System

clear all; 

N = 10^6; % number of messages for simulation % could also use 10^6 for more accurate results
n = 7; % codeword length
k = 4; % message length
R = k/n; % coding rate

% Generator Matrix %
P = [1 1 0; 0 1 1; 1 1 1; 1 0 1];
I_k = eye(k);
G = horzcat(P,I_k);

% Decoder Matrix %
I_n_k = eye(n-k);
H = horzcat(I_n_k,P');

% Decoding Error Table %
e1 = [0 0 0 0 0 0 0];
e2 = [0 0 1 0 0 0 0];
e3 = [0 1 0 0 0 0 0];
e4 = [0 0 0 0 1 0 0];
e5 = [1 0 0 0 0 0 0];
e6 = [0 0 0 0 0 0 1];
e7 = [0 0 0 1 0 0 0];
e8 = [0 0 0 0 0 1 0];
E = [e1;e2;e3;e4;e5;e6;e7;e8];

% Energy per bit definitions
Ec = 1; % energy in a coded bit
Eb = Ec/R; % energy in an uncoded bit

Ec_No_dB = [-3:1:10]; % range of Ec/N0 in dB
No = 1./(10.^(Ec_No_dB./10))*Ec; % Noise energy


for i = 1:length(Ec_No_dB) % iterating for various Ec/N0 values 

% Hamming Encoder %
u = randi([0,1],N,k); % generate random data
c = mod(u*G,2);  % generate Hamming codewords

% Modulator %
Tx = c; % copy values
Tx(Tx==0)=-1; % replace 0s with -1

% Channel (AWGN) %
% we can generate noise in one direction as this alone impacts demodulation
noise = sqrt(No(i)/2)*(randn(N,n)); % generate random noise
Rx = Tx + noise; % add to transmitted signal

% Reciever - Hard Decision %
r = Rx; % copy values
r(r<0)=0; % less than 0 assumed to be 0
r(r>0)=1; % greater than 0 assumed to be 1

% Hamming Decoder %
s = mod(r*H',2); % extract syndrome
s_index =  bin2dec(num2str(s))+1; % extract index of syndrome

% Correct Errors %
r_hat = mod(r + E(s_index,1:n),2); % correct errors by adding error vector
u_hat = r_hat(1:N,k:n); % extract estimated information bits

% Compute BER %
BER(i) = sum(u_hat(:) ~= u(:))/(N*k); % count errors and divide by total number of bits
BER_Theory(i) = qfunc(sqrt((2*Eb)/(No(i)))); % compute theoretical BER
end

% Plotting Results %
semilogy(10*log10(Eb./No),BER_Theory,'-b')
hold on
semilogy(10*log10(Eb./No),BER,'x-r')
ylim([10^-6 10^0])
xlim([0 10])
xlabel('E_b/N_0 (dB)')
legend('BPSK - Uncoded', 'BPSK - Hamming (7,4) Encoded ');
ylabel('BER')
grid on
hold off