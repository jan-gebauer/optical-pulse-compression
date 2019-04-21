%Working with 3.2.1 and 3.2.9

clc; clear all; close all;

b2 = -20;%-20e-27;
% alpha = -10000;%1e-5;
% gama = 0.01;%0.003;
To= 10;%125e-12; %Pulse width
T = 20*To;%8192e-12;
nSamples = 2^15;
v = linspace(0,nSamples-1,nSamples);
dT = T/nSamples;

Fs = (nSamples-1)/T;% dt = T/(nSamples-1)
t = -T/2:1/Fs:T/2;

df = 2*pi/T;

Fmax = Fs/2;
p = find(v > floor(nSamples/2));
v(p) = v(p)-nSamples;
f = v*df;
% f = 2*pi.*f;

dfx = 1/T;
fx = -Fmax:dfx:Fmax;
fx = 2*pi.*fx;
fx = fftshift(fx);

fy = (-nSamples/2:nSamples/2-1)*df;
fy = fftshift(fy);

f_vector(1,:) = f;
f_vector(2,:) = fx;
f_vector(3,:) = fy;

f = fy;

% Po = abs(b2)/(gama*To^2); %Peak power
Ld = To^2/abs(b2);
Lnl = Ld/10;
L = Ld;
%dz? %is this what Naks was talking about?
dz = min(Ld,Lnl)/100;
z_vector = 0:dz:2*L;

% Nz = 100;
% distance = 0.04;
% dz = distance/Nz;
% z_vector = 0:dz:distance-dz;

anal_wave = zeros(length(z_vector), length(t));
a = To/((To^2 - sqrt(-1)*b2*z_vector(1))^0.5);
b = exp(-(t.^2)/(2*(To^2 - sqrt(-1)*b2*z_vector(1))));
anal_wave(1,:) = abs(a*b).^2;
ToFWHM=find(abs(anal_wave(1,:))>abs(max(anal_wave(1,:))/2));
ToFWHM=length(ToFWHM);
b_factor(1) = ToFWHM/ToFWHM;
for i = 2:length(z_vector)
        a = To/((To^2 - sqrt(-1)*b2*z_vector(i))^0.5);
        b = exp(-(t.^2)/(2*(To^2 - sqrt(-1)*b2*z_vector(i))));
        anal_wave(i,:) = abs(a*b).^2;
        
        TFWHM=find(abs(anal_wave(i,:))>abs(max(anal_wave(i,:))/2));
        TFWHM=length(TFWHM);
        b_factor(i) = TFWHM/ToFWHM;
        
        

end
% b_factor = zeros(length(z_vector),1);
% ToFWHM=find(abs(anal_wave(1,:))>abs(max(anal_wave(1,:))/2));
% ToFWHM=length(ToFWHM);
% b_factor(1) = 1;
% for i = 2:length(z_vector)
%     TFWHM=find(abs(anal_wave(i,:))>abs(max(anal_wave(i,:))/2));
%     TFWHM=length(TFWHM);
%     b_factor(i) = TFWHM/ToFWHM;
% end
%%
simul_wave = zeros(length(z_vector), length(t));

C = 0;
chirp = 1+1i*C;
A = exp(-0.5*chirp*t.^2/To^2);
ToFWHM=find(abs(A)>abs(max(A)/2));
ToFWHM=length(ToFWHM);

D = -dz*(b2/2) * 1i * f.^2;

% f = ifftshift(f);
for i = 1:length(z_vector)
    At0 = fft(A);
    Atilde = At0.*exp(1i*b2*f.^2*dz/2); %At0.*exp(dz* (b2/2) * 1i * f.^2) 
    A = ifft(Atilde);
    simul_wave(i,:) = abs(A).^2;
    %b_factor(i) = To*sqrt([1+(z_vector(i)/Ld)^2]);
    
%     TFWHM=find(abs(A)>abs(max(A)/2));
%     TFWHM=length(TFWHM);
%     b_factor(i) = TFWHM/ToFWHM;
end

%%
% figure;
% mesh(t,z_vector, simul_wave); %(1:5) (1:5,:)
% xlabel('time');
% ylabel('distance');
% title('Simulated GVD-only gaussian in time domain');
% view([0 90]);
% 
% figure;
% mesh(t/To, z_vector, anal_wave); %(1:500,:)
% xlabel('time');
% ylabel('distance');
% title('Analytical GVD-only gaussian in time domain');

%%
figure;
plot(t/To,anal_wave(1,:));
%title('A');
%figure;
hold on;
plot(t/To,anal_wave(1000,:));
%title('B');
%figure;
hold on;
plot(t/To,anal_wave(end,:));
%title('C');
legend('z = 0', 'z = L/2', 'z = L')
% title('Analytical GVD-only Gaussian simulation in time domain')
xlabel('Time delay T/To')
ylabel('Intensity (a.u.)')

%%
figure;
plot(t/To,simul_wave(1,:));
%title('A');
%figure;
hold on;
plot(t/To,simul_wave(1000,:));
%title('B');
%figure;
hold on;
plot(t/To,simul_wave(end,:));
%title('C');
legend('z = 0', 'z = L/2', 'z = L')
title('Numerical GVD-only Gaussian simulation in time domain')
xlabel('Time delay T/To')
ylabel('Amplitude')

%%
Xf_ind = t(find(simul_wave(1,:) == max(simul_wave(1,:))));
Xf = max(simul_wave(1,:));
Xth_ind = t(find(simul_wave(1000,:) == max(simul_wave(1000,:))));
Xth = max(simul_wave(1000,:));
%%
figure;
plot(z_vector,b_factor)
% title('Broadening factor for GVD-only Gaussian')
ylabel('Broadening factor')
xlabel('Distance z (km)')