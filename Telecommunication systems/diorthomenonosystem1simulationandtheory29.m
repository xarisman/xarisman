N = 10^5; %number of symbols to process
M = 4;
a = 0.25;
alpha4pam = [-3*sqrt(a/5) -sqrt(a/5) sqrt(a/5) 3*sqrt(a/5)]; %cordinates of 4-PAM
alpha4pam2 = [-3*sqrt((1-a)/5) -sqrt((1-a)/5) sqrt((1-a)/5) 3*sqrt((1-a)/5)];
SNR = -3:20;
ipHat = zeros(1,N); %creates array with N datatype

%Simulation for U1
for ii = 1:length(SNR)
    ip = randsrc(1, N, alpha4pam);
    s = (1/sqrt(a))*ip; %normaliazation of energy
    n = 1/sqrt(a*4)*(randn(1,N) + (1j)*randn(1,N));
    
    y = s + 10^(-SNR(ii)/20)*n; % additive white gaussian noise

    %demodulation
    r = real(y);
    %areas of decisions
    ipHat(find(r< -2/sqrt(5))) = -3*sqrt(a/5);
    ipHat(find(r>= 2/sqrt(5))) = 3*sqrt(a/5);
    ipHat(find(r>=-2/sqrt(5) & r<0)) = -sqrt(a/5);
    ipHat(find(r>=0 & r<2/sqrt(5))) = sqrt(a/5);
    
    nErr(ii) = size(find([ip - ipHat]), 2); %counting errors
end

%Simulation for U2
for ii = 1:length(SNR)
    ip = randsrc(1, N, alpha4pam2);
    s = (1/sqrt(1-a))*ip; %normaliazation of energy
    n = 1/sqrt((1-a)*8)*(randn(1,N) + (1j)*randn(1,N));
    
    y = s + 10^(-SNR(ii)/20)*n; % additive white gaussian noise

    %demodulation
    r = real(y);
    %areas of decisions
    ipHat(find(r< -2/sqrt(5))) = -3*sqrt((1-a)/5);
    ipHat(find(r>= 2/sqrt(5))) = 3*sqrt((1-a)/5);
    ipHat(find(r>=-2/sqrt(5) & r<0)) = -sqrt((1-a)/5);
    ipHat(find(r>=0 & r<2/sqrt(5))) = sqrt((1-a)/5);
    
    nErr2(ii) = size(find([ip - ipHat]), 2); %counting errors
end

simBer = nErr/N; %propabulity of error;
simBer2 = nErr2/N;
%Theoretical results for the 2 users
theorySer = 0.75*erfc(sqrt(0.2*a*2*(10.^(SNR/10)))); %forgot a
theorySer2 = 0.75*erfc(sqrt(0.2*4*(1-a)*(10.^(SNR/10))));
close all
figure
semilogy(SNR, theorySer, 'b.-');
hold on
semilogy(SNR, simBer, 'mx-');
hold on
semilogy(SNR, theorySer2, 'g.-');
hold on
semilogy(SNR, simBer2, 'k.-');

grid on
legend('theory', 'simulation');
xlabel('SNR, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for 4-PAM modulation')