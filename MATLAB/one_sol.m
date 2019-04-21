%First order soliton with chirp
%%
clc; clear all; close all;

b2 = -20;% -10;
b2_initial = b2;
gama = abs(b2)/100;

To= 10; %Pulse width
T = 20*To;
nSamples = 2^3;

% C = 0.5;
C= 0;

Fs = (nSamples-1)/T;

dt = 1/Fs;
t = (-nSamples/2:nSamples/2-1)*dt;

df = 2*pi/T;
f = (-nSamples/2:nSamples/2-1)*df;
f = fftshift(f);

N_order = 1;
Po = N_order^2;
Ld = To^2/abs(b2);
Lnl = 1/(Po*gama);
L = 50;


dz = min(Ld,Lnl)/100;
z_vector = 0:dz:L;

%%
simul_wave = zeros(length(z_vector), length(t));
spec_wave = zeros(length(z_vector), length(t));
A = N_order*sech(t/To).*exp(1i*C*t.^2/(2*To^2));
% A = exp(-0.5*(1+1i*C)*t.^2/To^2);
fwhm1=find(abs(A)>abs(max(A)/2));
fwhm1=length(fwhm1);
fwhm_vec = zeros(1, length(z_vector));
fwhm_vec(1) = To;
%%
% f = ifftshift(f);
b2 = b2_initial;
for i = 1:length(z_vector)
    D = exp(1i*(dz/2)*b2*f.^2).*fft(A);
    D = ifft(D);
    N = exp(1i*gama*abs(A).^2*dz);
    A = D.*N;
%     spec_wave(i,:) = ifftshift(abs(fft(A)));
%     test_wave(i,:) = A;
    simul_wave(i,:) = abs(A).^2;
    

%     plot(t/To, simul_wave(i,:));
%     ylim([0 1.5])
%     pause(.00000001);
    
    fwhm=find(abs(A)>abs(max(A)/2));
    fwhm=length(fwhm);
    fwhm_vec(i) = fwhm;
    b_factor(i) = fwhm/fwhm1;
end

%%
% figure;
% mesh(t/To,z_vector, simul_wave);
% xlabel('Time Delay T/To');
% ylabel('Distance z (km)');
% zlabel('Intensity (a.u.)');
% % title('First order soliton with initial chirp');
% view(0, 90);
%%
figure;
plot(t/To,simul_wave(1,:));
hold on;
% plot(t/To, simul_wave(2,:));
% hold on;
plot(t/To, simul_wave(end,:), 'o');
legend('z = 0', 'z = L');
% title({'Soliton pulse', 'in the anomalous dispersion regime'})
xlabel('Time Delay T/To')
ylabel('Intensity (a.u.)')

%%
% f = ifftshift(f);
% %%
% figure;
% plot(f,spec_wave(1,:));
% hold on;
% plot(f,spec_wave(80,:));
% hold on;
% plot(f,spec_wave(end,:));
% legend('z = 0', 'z = L/2', 'z = L');
% title('Frequency spectrum')
% xlabel('f')
% ylabel('Amplitude')
% xlim([-5 5])
%%
figure;
plot(z_vector,b_factor)
% title('Broadening factor')
xlabel('Distance z (km)');
ylabel('Broadening factor');
%%
% eT = trapz(simul_wave(end,:));
% simAmp = max(simul_wave(end,:));
% finWidth = fwhm_vec(end);
% eSech = 2*simAmp*(finWidth/1.763); %1.763 is a constant from Cao Wai 2005 paper
% PE = abs(eT-eSech)/eT *100;
% checkN = (gama*simAmp*To^2)/abs(b2);
%%
% figure
% subplot(1,2,1)
% mesh(t/To,z_vector, simul_wave);
% xlabel('Time Delay T/To');
% ylabel('Distance z');
% zlabel('Amplitude');
% view([0 90])
% subplot(1,2,2)
% plot(t/To,simul_wave(1,:));
% hold on;
% plot(t/To, simul_wave(end,:), 'o');
% xlabel('Time Delat T/To');
% ylabel('Amplitude');
% legend('z = 0', 'z = L');
%%
% figure
% subplot(1,3,1)
% mesh(t/To,z_vector, simul_wave);
% xlabel('Time Delay T/To');
% ylabel('Distance z');
% zlabel('Amplitude');
% view([0 90])
% subplot(1,3,2)
% plot(t/To,simul_wave(1,:));
% hold on;
% plot(t/To, simul_wave(end,:), 'o');
% %legend('z = 0', 'z = L');
% xlabel('Time Delay T/To');
% ylabel('Amplitude');
% subplot(1,3,3)
% plot(z_vector,b_factor)
% xlabel('Distance z');
% ylabel('Broadening Factor');
%%
% c_imp = readmatrix('vals.csv');
% c_imp = c_imp(:,1:end-1);
%%
% figure
% mesh(c_imp)
% %%
% figure
% plot(abs(A).^2);
% hold on;
% plot(c_imp(1,:))
% legend('A','C');