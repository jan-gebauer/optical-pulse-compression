%First order soliton with exponential DDF
%%
clc; clear all;

b2 = -33;%-20; %GVD
b2_initial = b2;

sigma = 0.05; %Decay rate for the exponential function
% sigma = C*b2_initial;
% C = sigma/b2;
C= 0;

%Time, and frequency vectors
To= 10; %Pulse width
T = 20*To;
nSamples = 2^15;
Fs = (nSamples-1)/T;

dt = 1/Fs;
t = (-nSamples/2:nSamples/2-1)*dt;

df = 2*pi/T;
f = (-nSamples/2:nSamples/2-1)*df;
f = fftshift(f);

%Additional setup
N_order = 1;
Po = N_order^2;
Ld = To^2/abs(b2);
Lnl = Ld;
gama = 1/(Po*Lnl);
L = 30;%10*Ld;

%distance vectors
dz = min(Ld,Lnl)/100;
z_vector = 0:dz:L;

%%
%Initial Wave setup
simul_wave = zeros(length(z_vector), length(t));
spec_wave = zeros(length(z_vector), length(t));
A = N_order*sech(t/To).*exp(1i*C*t.^2/(2*To^2));
Ao = abs(A).^2;
%%
fwhm1=find(abs(Ao)>abs(max(Ao)/2));
fwhm1=length(fwhm1);
fwhm_vec = zeros(1, length(z_vector));
%%
%Split Step Fourier Method
b2 = b2_initial;

b2_vec = b2_initial;
for i = 1:length(z_vector)
    D = exp(1i*(dz/2)*b2*f.^2).*fft(A);
    D = ifft(D);
    N = exp(1i*gama*abs(A).^2*dz);
    A = D.*N;
    spec_wave(i,:) = ifftshift(abs(fft(A)));
    simul_wave(i,:) = abs(A).^2;
    

%     plot(t/To, simul_wave(i,:));
%     ylim([0 4])
%     pause(.00000001);
   
    
    b2_vec(i) = b2;
    b2 = b2_initial*exp(-sigma*z_vector(i));
    
    fwhm=find(abs(simul_wave(i,:))>abs(max(simul_wave(i,:))/2));
    fwhm=length(fwhm);
    fwhm_vec(i) = fwhm;
    c_factor(i) = fwhm1/fwhm;
end

%%
%This is the 3D plot
% figure;
% mesh(t/To,z_vector, simul_wave);
% xlabel('Time Delay T/To');
% ylabel('Distance z');
% zlabel('Amplitude');
% title('First order soliton with initial chirp');
% view(0, 90);
%%
%2D Plot
figure;
plot(t/To,simul_wave(1,:));
hold on;
plot(t/To, simul_wave(end,:));
% legend('z = 0', 'z = L');
% title({'Soliton pulse', 'in the anomalous dispersion regime'})
xlabel('Time Delay T/To')
ylabel('Intensiry (a.u.)')

%%
% f = ifftshift(f);
%%
% figure;
% plot(f,spec_wave(1,:));
% hold on;
% plot(f,spec_wave(20,:));
% hold on;
% plot(f,spec_wave(end,:));
% legend('z = 0', 'z = L/2', 'z = L');
% title('Frequency spectrum')
% xlabel('f')
% ylabel('Amplitude')
% xlim([-5 5])

% %%
% figure
% plot(b2_vec)
%%
% au_wave = simul_wave(1,:);%/min(simul_wave(end,:));
% au_wave = normalize(au_wave,'DataVariables',simul_wave(end,:));
%plot(log(au_wave)/log(2))

%%
for i = 1:length(z_vector)
    eT(i) = trapz(t/To,simul_wave(i,:));
    simAmp(i) = max(simul_wave(i,:));
    eSech(i) = 2*simAmp(i)*((fwhm_vec(i)*dt/To)/1.763); %1.763 is a constant from Cao Wai 2005 paper
    PE(i) = abs(eT(i)-eSech(i))/eT(i) *100;
end
%%
figure;
plot(z_vector,c_factor)
% title('Compression Factor')
xlabel('Distance z (km)');
ylabel('Compression Factor')
% hold on
% plot(z_vector,PE)
%%
fake_z_vector = 1:L;
for i = 1:L
    tar = find(z_vector < fake_z_vector(i)+dz/2 & z_vector > fake_z_vector(i)-dz/2);
    fake_sim_wave(i,:) = simul_wave(tar,:);    
end
%%
figure
surface(t/To, fake_z_vector, fake_sim_wave, 'FaceAlpha', [0], 'MeshStyle', 'row');
xlim([-2 2])
grid on
grid minor
view([-45 30])
xlabel('Time Delay T/To');
ylabel('Distance z');
zlabel('Amplitude');