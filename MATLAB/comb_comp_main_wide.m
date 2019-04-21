%First order soliton with chirp
%%
clc; clear; close all;

b2 = -200;%-33;% -10;
b2_initial = b2;

To= 10;%0.908;%10; %Pulse width
T = 4*20*To;
nSamples = 2^15;

sigma = 5.3;%0.05;
C = sigma/b2;

% C = -0.265;
% sigma = C*b2;

Fs = (nSamples-1)/T;

dt = 1/Fs;
t = (-nSamples/2:nSamples/2-1)*dt;

df = 2*pi/T;
f = (-nSamples/2:nSamples/2-1)*df;
f = fftshift(f);

N_order = 1;
Ld = To^2/abs(b2);
gama = 20;
L = 0.8;
Po = 0.25*abs(b2)/gama/To^2;%3;
Lnl = 1/(Po*gama);


dz = min(Ld,Lnl)/100;
z_vector = 0:dz:L;

sep_factor = 0.5*To;
%%
simul_wave = zeros(length(z_vector), length(t));
% spec_wave = zeros(length(z_vector), length(t));

%Three pulses
% A = sech(t/To + -1*sep_factor*To).*exp((1i*C*t.^2)/2);
% A = A + sech(t/To + 0*sep_factor*To).*exp((1i*C*t.^2)/2);
% A = A + sech(t/To + 1*sep_factor*To).*exp((1i*C*t.^2)/2);

%Two pulses
% A = sech(t/To + sep_factor/2 + -1*sep_factor).*exp((1i*C*t.^2)/2);
% A = A + sech(t/To + sep_factor/2 + 0*sep_factor).*exp((1i*C*t.^2)/2);

%Four pulses
A = sech(t/To + sep_factor/2 + -2*sep_factor).*exp((1i*C*t.^2)/2);
A = A + sech(t/To + sep_factor/2 + -1*sep_factor).*exp((1i*C*t.^2)/2);
A = A + sech(t/To + sep_factor/2 + 0*sep_factor).*exp((1i*C*t.^2)/2);
A = A + sech(t/To + sep_factor/2 + 1*sep_factor).*exp((1i*C*t.^2)/2);


Ao = abs(A).^2;
plot(t/To,Ao)
fwhm1=find(abs(Ao)>abs(max(Ao)/2));
fwhm1=length(fwhm1);
fwhm_vec = zeros(1, length(z_vector));

eT = 0;
sinAmp = 0;
finWidth = 0;
eSech = 0; %1.763 is a constant from Cao Wai 2005 paper
PE = 0;
checkN = 0;

%%
% figure;
% mesh(t/To,z_vector, simul_wave);
% % mesh(t/To,z_vector(1:511),simul_wave(1:511,:))
% xlabel('Time Delay T/To');
% ylabel('Distance z/LD');
% zlabel('Amplitude');
% title('First order soliton with initial chirp');

%%
% A = A_vec(511,:);
% Ao = abs(A).^2;
% % plot(t/To,Ao)
% fwhm1=find(abs(Ao)>abs(max(Ao)/2));
% fwhm1=length(fwhm1);
% fwhm_vec = zeros(1, length(z_vector));
%%
for i = 1:length(z_vector)

    D = exp(1i*(dz/2)*b2*f.^2).*fft(A);
    D = ifft(D);
    N = exp(1i*gama*abs(A).^2*dz);
    A = D.*N;
%     spec_wave(i,:) = ifftshift(abs(fft(A)));
    simul_wave(i,:) = abs(A).^2;
%     plot(t/To, simul_wave(i,:));
%     ylim([0 2])
%     pause(.00000001);
    
%     fwhm=find(abs(A)>abs(max(A)/2));
%     fwhm=length(fwhm);
% %     fwhm_vec(i) = fwhm;
%     fwhm_vec = [fwhm_vec fwhm];
%     c_factor(i) = fwhm/fwhm1;
    
    
    fwhm=find(abs(A).^2>max(abs(A).^2)/2);
    fwhm=length(fwhm);
    fwhm_vec(i) = fwhm;
    c_factor(i) = (fwhm1/fwhm);
    
    b2_vec(i) = b2;
    b2 = b2_initial*exp(-sigma*z_vector(i));

    
    
    eT = [eT trapz(simul_wave(i,:))];
    simAmp = [sinAmp max(simul_wave(i,:))];
    finWidth = [finWidth fwhm_vec(i)];
    eSech = [eSech 2*simAmp(end)*(finWidth(end)/1.763)]; %1.763 is a constant from Cao Wai 2005 paper
    PE = [PE abs(eT(end)-eSech(end))/eT(end) *100];
    checkN = [checkN (gama*simAmp(end)*To^2)/abs(b2)];
end

%%
figure;
mesh(t/To,z_vector, simul_wave);
xlabel('Time Delay T/To');
ylabel('Distance z/LD');
zlabel('Amplitude');
title('First order soliton with initial chirp');
% view(0, 90);
% %%
% figure;
% plot(t/To,simul_wave(1,:));
% hold on;
% plot(t/To, simul_wave(2,:));
% hold on;
% plot(t/To, simul_wave(end,:));
% legend('z = 0', 'z = L/2', 'z = L');
% title('Soliton pulses')
% xlabel('T/To')
% ylabel('Amplitude')

%%
figure;
plot(PE)
hold on
plot(c_factor)
xlabel('Distance z');
ylabel('Compression Factor');
%%
% eT = trapz(simul_wave(end,:));
% simAmp = max(simul_wave(end,:));
% finWidth = fwhm_vec(end);
% eSech = 2*simAmp*(finWidth/1.763); %1.763 is a constant from Cao Wai 2005 paper
% PE = abs(eT-eSech)/eT *100;
% checkN = (gama*simAmp*To^2)/abs(b2);
%%
figure
plot(z_vector,b2_vec)

%%
fake_z_vector = 1/100:1/100:L;
for i = 1:(L*100-1)
   pos = find(z_vector > fake_z_vector(i)-(dz/2), 1);
   fake_sim_wave(i,:) = simul_wave(pos,:);    
end
fake_sim_wave(8,:) = simul_wave(end,:);
%%
figure
surface(t/To, fake_z_vector, fake_sim_wave, 'FaceAlpha', [0], 'MeshStyle', 'row');
xlim([-6 6])
grid on
grid minor
view([-45 30])
xlabel('Time Delay T/To');
ylabel('Distance z');
zlabel('Amplitude');