clear all;
N = 10^5; %number of symbols to process
th = 1:44;
%theta = 26;
p = 0.5;
SNR = 1:20;
ipHat = zeros(1,N);
ipHat2 = zeros(1,N);%creates array with N datatype
d = sqrt(p);
d2 = sqrt(1-p);
array_3 = zeros(1,length(th));
array_2 = zeros(1,length(SNR));
array_1 = zeros(1,length(SNR));

for i = 1:length(SNR)
    for theta = 1:length(th)
        alpha4pam = [-cosd(45-theta)*d -sind(45-theta)*d sind(45-theta)*d cosd(45-theta)*d];
        alpha4pam2 = [-cosd(45-theta)*d2 -sind(45-theta)*d2 sind(45-theta)*d2 cosd(45-theta)*d2];
        
        %Simulation of U1
        ip = randsrc(1, N, alpha4pam);
        s = ip/sqrt(1/(2*p));
        n = 1/sqrt(32*p)*(randn(1,N) + (1j)*randn(1,N));
        
        y = s + 10^(-i/20)*n;
    
        r = real(y);
        ipHat(find(r< -(cosd(45-theta)*d+sind(45-theta)*d)/2*sqrt(2*p))) = -cosd(45-theta)*d;
        ipHat(find(r>= (cosd(45-theta)*d+sind(45-theta)*d)/2*sqrt(2*p))) = cosd(45-theta)*d;
        ipHat(find(r>=-(cosd(45-theta)*d+sind(45-theta)*d)/2*sqrt(2*p) & r<0)) = -sind(45-theta)*d;
        ipHat(find(r>=0 & r<(cosd(45-theta)*d+sind(45-theta)*d)/2*sqrt(2*p))) = sind(45-theta)*d;
        nErr = size(find([ip - ipHat]), 2); %counting errors
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Simulation of U2
        ip2 = randsrc(1, N, alpha4pam2);
        s2 = ip2./sqrt(1/(2*(1-p)));
        n2 = 1/sqrt(64*(1-p))*(randn(1,N) + (1j)*randn(1,N));
        
        y2 = s2 + 10^(-i/20)*n2;
    
        r2 = real(y2);
        ipHat2(find(r2< -(cosd(45-theta)*d2+sind(45-theta)*d2)/2*sqrt(2*(1-p)))) = -cosd(45-theta)*d2;
        ipHat2(find(r2>= (cosd(45-theta)*d2+sind(45-theta)*d2)/2*sqrt(2*(1-p)))) = cosd(45-theta)*d2;
        ipHat2(find(r2>=-(cosd(45-theta)*d2+sind(45-theta)*d2)/2*sqrt(2*(1-p)) & r2<0)) = -sind(45-theta)*d2;
        ipHat2(find(r2>=0 & r2<(cosd(45-theta)*d2+sind(45-theta)*d2)/2*sqrt(2*(1-p)))) = sind(45-theta)*d2;
        nErr2 = size(find([ip2 - ipHat2]), 2);
           
        array_3(theta) = max(nErr/N, nErr2/N);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Ser2(i) = 1/2*erfc(sqrt((sind(th).^2).*8*(1-p)*(10.^(i/10))))+1/4.*erfc(sqrt((cosd(th)-sind(th)).^2.*8*(1-p)*(10.^(i/10))));
    %Ser1(i) = 1/2*erfc(sqrt((sind(th).^2).*4*p*(10.^(i/10))))+1/4.*erfc(sqrt((cosd(th)-sind(th)).^2.*4*p*(10.^(i/10))));
    end

    %Saving the optimazid theta
    semilogy(th, array_3)
    hold on
    [M, I] = min(array_3);
    array_1(i) = I;
    
end

%{
simBer = nErr/N;
simBer2 = nErr2/N;
semilogy(simBer, 'g')
hold on
semilogy(simBer2, 'y')
hold on
plot(Ser2, 'b')
hold on
plot(Ser1, 'r')
%}