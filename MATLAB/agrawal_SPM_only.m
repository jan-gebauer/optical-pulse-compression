%Working with 4.1.1 and 4.1.4
%%
clc; clear all; close all;

b2 = 50;%-20e-27;
alpha = -1;%1e-5;
gama = 2;%0.003;
To= 1;%125e-12; %Pulse width
T = 20*To;%8192e-12;
nSamples = 2^10;
%%

Fs = (nSamples-1)/T;% dt = T/(nSamples-1)
Fmax = Fs/2;
df = 1/T;
t = -T/2:1/Fs:T/2;
f = -Fmax:df:Fmax;
f = 2*pi.*f;
%%
% Po = abs(b2)/(gama*To^2); %Peak power
Po = 1;
Ld = 10000*To^2/abs(b2);
Lnl = (gama*Po)^-1;
L = 10*Lnl; %There is no point in doing this, is there? As crazy as it sounds
%dz? %is this what Naks was talking about?
dz = min(Ld,Lnl)/100;
z_vector = 0:dz:L;
%%

anal_wave = zeros(length(z_vector), length(t));
anal_wave_time = zeros(length(z_vector), length(t));
% % U = U*exp(1i*phi); %U is not a constant at all
Uo = exp(-0.5*t.^2/To^2);
% % Leff = (1 - exp(-alpha*L))/alpha; %is this a constant? What is L?
% % phi = abs(Uo).^2 * Leff/Lnl; %phi is a constant in terms of distance
%%
for i = 1:length(z_vector)
    if alpha == 0
        Leff = L;
    else
        Leff = (1 - exp(-alpha*z_vector(i)))/alpha;
    end
    phi = abs(Uo).^2 .* Leff/Lnl;
    U = Uo.*exp(1i*phi);
    U_fft = fftshift(fft(U))/Fs;
    anal_wave(i,:) = abs(U_fft).^2;
    anal_wave_time(i,:) = abs(U).^2;
end
%%
simul_wave = zeros(length(z_vector), length(t));
simul_wave_time = zeros(length(z_vector), length(t));
A = exp(-0.5*t.^2/To^2);
Ao = A;
A_fft = fftshift(fft(A))/Fs;
simul_wave(1,:) = abs(A_fft).^2;
simul_wave_time(1,:) = abs(A).^2;
%%
f = ifftshift(f);
for i = 2:length(z_vector)
    N = (1i*exp(-alpha*z_vector(i)) / Lnl)*abs(A).^2;
    A = exp(dz*N).*A;
    A_fft = fftshift(fft(A))/Fs;
    simul_wave(i,:) = abs(A_fft).^2;
    simul_wave_time(i,:) = abs(A).^2;
end
f = fftshift(f);
%%
% figure;
mesh(f,z_vector, simul_wave); %(1:5) (1:5,:)
xlabel('time');
ylabel('distance');
title('Simulated GVD-only gaussian in time domain');

figure;
mesh(f, z_vector, anal_wave); %(1:500,:)
xlabel('time');
ylabel('distance');
title('Analytical GVD-only gaussian in time domain');
%%
% close all;
figure;
subplot(3,2,2);
plot(t/To,anal_wave_time(1,:));
title('z = 0')
subplot(3,2,4);
plot(t/To,anal_wave_time(400,:));
title('z = L/2')
ylabel('Intensity (a.u.)')
subplot(3,2,6);
plot(t/To,anal_wave_time(533,:));
title('z = 2L/3')
xlabel('Time Delay T/To')

subplot(3,2,1);
plot(f,anal_wave(1,:));
title('z = 0')
xlim([-100 100]);
subplot(3,2,3);
plot(f,anal_wave(400,:));
title('z = L/2')
ylabel('Intensity (a.u.)')
xlim([-100 100]);
subplot(3,2,5);
plot(f,anal_wave(533,:));
title('z = 2L/3')
xlabel('Frequency \omega/Hz')
xlim([-100 100]);
%%
figure;
subplot(3,2,2);
plot(t/To,simul_wave_time(1,:));
title('z = 0')
subplot(3,2,4);
plot(t/To,simul_wave_time(400,:));
title('z = L/2')
ylabel('Amplitude')
subplot(3,2,6);
plot(t/To,simul_wave_time(533,:));
title('z = 2L/3')
xlabel('Time Delay T/To')

subplot(3,2,1);
plot(f,simul_wave(1,:));
title('z = 0')
xlim([-100 100]);
subplot(3,2,3);
plot(f,simul_wave(400,:));
title('z = L/2')
ylabel('Amplitude')
xlim([-100 100]);
subplot(3,2,5);
plot(f,simul_wave(533,:));
title('z = 2L/3')
xlabel('Frequency \omega/Hz')
xlim([-100 100]);