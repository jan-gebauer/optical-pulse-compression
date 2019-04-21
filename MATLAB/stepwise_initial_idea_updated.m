%Stepwise approximation of exponential DDF
%Add Sigma
%Change the stepping

%%
clc; clear all; close all;

b2 = -33;% -10;
b2_initial = b2;

To= 10; %Pulse width
T = 20*To;
nSamples = 2^13;

sigma = 0.05;
C = sigma/b2;

Fs = (nSamples-1)/T;

dt = 1/Fs;
t = (-nSamples/2:nSamples/2-1)*dt;

df = 2*pi/T;
f = (-nSamples/2:nSamples/2-1)*df;
f = fftshift(f);

N_order = 1;
Po = N_order^2;
Ld = To^2/abs(b2);
Lnl = Ld;
gama = 1/(Po*Ld);
L = 30;

dz = min(Ld,Lnl)/100;
z_vector = 0:dz:L;

%%
simul_wave = zeros(length(z_vector), length(t));
spec_wave = zeros(length(z_vector), length(t));
A = N_order*sech(t/To).*exp(1i*C*t.^2/(2*To^2));
fwhm1=find(abs(abs(A).^2)>abs(max(abs(A).^2)/2));
fwhm1=length(fwhm1);
fwhm_vec = zeros(1, length(z_vector));
%%
% f = ifftshift(f);

tot_steps = length(z_vector); %why 63?
segs = 4; %input variable
b2_step = tot_steps/segs;
b2_vec_inc = 0;
b2_vec = b2_initial;
for i = 1:length(z_vector)

%     b2 = b2_initial+((-b2_initial + b2_final)/(length(z_vector)))*i;
%     if i >= taper_start + taper_len + 1
%         gama = -b2/100;
%     end
    D = exp(1i*(dz/2)*b2*f.^2).*fft(A);
    D = ifft(D);
    N = exp(1i*gama*abs(A).^2*dz);
    A = D.*N;
    spec_wave(i,:) = ifftshift(abs(fft(A)));
    simul_wave(i,:) = abs(A).^2;
    

%     plot(t/To, simul_wave(i,:));
%     ylim([0 1.5])
%     pause(.00000001);

    
%     if i >= taper_start && i < taper_start + taper_len
%         b2_vec(i) = b2;
%         b2 = b2_initial+((-b2_initial + b2_final)/(taper_len))*i;
%         b2 = b2_initial*exp(-0.032*(i - taper_start));
%         b2 = b2_initial+((-b2_initial + b2_final)/(length(z_vector)))*i;
%         b2
        
        if fix(mod((i),tot_steps/segs)) == 0
%             b2_vec_inc = b2_vec_inc + b2_step;
%             b2_vec = [b2_vec b2_initial*exp(-0.032*b2_vec_inc)]; %???
            b2 = b2_initial*exp(-sigma*z_vector(i)); 
            b2_vec = [b2_vec b2];
%             b2 = b2_initial*exp(-0.032*b2_vec_inc);
        else
            b2_vec = [b2_vec b2];
        end

%         gama = (N_order^2 * abs(b2))/(max(simul_wave(i,:))*fwhm_vec(i)^2);
%         gama_vec = [gama_vec gama];
%         n_vec = [n_vec (gama*max(simul_wave(i,:))*fwhm_vec(i)^2)/abs(b2)];
%     end
    b2_vec_cont(i) = b2_initial*exp(-sigma*z_vector(i));
        
    fwhm=find(abs(abs(A).^2)>abs(max(abs(A).^2)/2));
    fwhm=length(fwhm);
    fwhm_vec(i) = fwhm;
    c_factor(i) = (fwhm1/fwhm);
    
%     if i ~= 1 && c_factor(i) < c_factor(i-1)
%         break
%     end

end

%%
figure;
mesh(t/To,z_vector, simul_wave);
xlabel('Time Delay T/To');
ylabel('Distance z/LD');
zlabel('Amplitude');
title('First order soliton with initial chirp');
% view(0, 90);
%%
% figure;
% plot(t/To,simul_wave(1,:));
% hold on;
% plot(t/To, simul_wave(20,:));
% hold on;
% plot(t/To, simul_wave(end,:));
% legend('z = 0', 'z = L/2', 'z = L');
% title({'Soliton pulse', 'in the anomalous dispersion regime'})
% xlabel('T/To')
% ylabel('Amplitude')
% 
% %%
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
hold on
plot(z_vector,PE)
title('compression factor')

% max_cf = max(c_factor);
%%
% eT = trapz(t/To,simul_wave(end,:));
% simAmp = max(simul_wave(end,:));
% finWidth = fwhm_vec(end);
% eSech = 2*simAmp*((finWidth*dt/To)/1.763); %1.763 is a constant from Cao Wai 2005 paper
% PE = abs(eT-eSech)/eT *100;

%%
figure
plot(b2_vec(1:end-1))
hold on;
plot(b2_vec_cont)
% title('Initial idea'); %how to actually call this?
xlabel('Distance z');
ylabel('Dispersion coefficient \beta_2');
% xlim([0 30])
% hold on;
% i_vec = 1:length(z_vector);
% i_vec_2 = [i_vec, fliplr(i_vec)];
% inBetween = [b2_vec(1:end-1), fliplr(b2_vec_cont)];
% fill(i_vec_2, inBetween, 'g');

%%
area_diff = trapz(abs(b2_vec)) - trapz(abs(b2_vec_cont));