%First order soliton with chirp
%%
clc; clear all; close all;

b2 = -20;% -10;
b2_initial = b2;

To= 10; %Pulse width
T = 20*To;
nSamples = 2^3;

C = -0.01;
ao = -4/((To^2) * pi);
% sigma = C*b2_initial;
sigma = ao*b2_initial;

%sigma = 0.05;
%C = sigma/b2;

Fs = (nSamples-1)/T;

dt = 1/Fs;
t = (-nSamples/2:nSamples/2-1)*dt;

df = 2*pi/T;
f = (-nSamples/2:nSamples/2-1)*df;
f = fftshift(f);

N_order = 2;
Po = 0.4;
Ld = To^2/abs(b2);
% Lnl = 1/(Po*gama);
Lnl = Ld;
gama = 2;
L = 12;%2*Ld;

dz = min(Ld,Lnl)/2000;
z_vector = 0:dz:L;
%%
b2_vec_cont = abs(b2_initial*exp(-sigma*z_vector));

vec = simpsonFunc(b2_vec_cont, z_vector,0.029);
length(vec)

for i = 1:length(vec)
    vec(i) = find(z_vector == vec(i), 1);
end
vec = [0 vec length(z_vector)];

seg_vec2 = [];
for i = 1:length(vec)-1
   segment_len =  vec(i+1) - vec(i);
   segment = ones(1,fix(segment_len))*b2_vec_cont(fix(segment_len/2 + vec(i)));
   seg_vec2 = [seg_vec2 segment];
end
seg_vec2 = [seg_vec2 ones(1,length(z_vector)-length(seg_vec2))*seg_vec2(end)];

%%
figure
plot(z_vector,(-1)*seg_vec2(1:length(z_vector)));
hold on
plot(z_vector,(-1)*b2_vec_cont);
title('Automated');
xlabel('Distance z');
ylabel('Dispersion coefficient \beta_2');
seg_vec2 = -1*seg_vec2;
%%
simul_wave = zeros(length(z_vector), length(t));
spec_wave = zeros(length(z_vector), length(t));
% A = N_order*sech(t/To).*exp(1i*C*t.^2/(2*To^2));
A = (2/To) * sqrt(b2_initial/gama) * sech(t/To) .* exp(1i*((ao*t.^2)/2 + (pi/8)));
% A = 2 * sech(t/To) .* exp(1i*((ao*t.^2)/2 + (pi/8)));
fwhm1=find(abs(abs(A).^2)>abs(max(abs(A).^2)/2));
fwhm1=length(fwhm1);
fwhm_vec = zeros(1, length(z_vector));
%%
% f = ifftshift(f);
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
%     ylim([0 30])
%     pause(.00000001);
   
    

    b2_vec(i) = b2;
%     b2 = b2_initial*exp(-sigma*z_vector(i));  
    b2 = seg_vec2(i);
    
    fwhm=find(abs(abs(A).^2)>abs(max(abs(A).^2)/2));
    fwhm=length(fwhm);
    fwhm_vec(i) = fwhm;
    c_factor(i) = (fwhm1/fwhm);
    
end

%%
% figure;
% mesh(t/To,z_vector, simul_wave);
% xlabel('Time Delay T/To');
% ylabel('Distance z');
% zlabel('Amplitude');
% title('First order soliton with initial chirp');
% view(0, 90);
%%
figure;
plot(t/To,simul_wave(1,:));
hold on;
plot(t/To, simul_wave(end,:));
xlabel('Time Delay T/To')
ylabel('Intensity (a.u.)')




%%
for i = 1:length(z_vector)
    eT(i) = trapz(t/To,simul_wave(i,:));
    simAmp(i) = max(simul_wave(i,:));
    eSech(i) = 2*simAmp(i)*((fwhm_vec(i)*dt/To)/1.763); %1.763 is a constant from Cao Wai 2005 paper
    PE(i) = abs(eT(i)-eSech(i))/eT(i) *100;
end
%%
fake_z_vector = [1:L];
for i = 1:L
   fake_sim_wave(i,:) = simul_wave(find(z_vector == fake_z_vector(i)),:);    
end
%%
figure
surface(t/To, fake_z_vector, fake_sim_wave, 'FaceAlpha', [0], 'MeshStyle', 'row');
xlim([-2 2])
grid on
grid minor
view([-45 30])
xlabel('Time Delay T/To');
ylabel('Distance z (km)');
zlabel('Intensity (a.u.)');
%%
figure;
plot(z_vector,c_factor)
xlabel('Distance z (km)');
ylabel('Compression Factor');
% title('breather compression factor')
hold on
plot(z_vector,PE)
% plot(exp(sigma*z_vector))


%%
function end_vec = simpsonFunc(exp_vec, distance_vec,tar)

%     figure
%     plot(distance_vec,exp_vec)
%     trapz(distance_vec,exp_vec);
    end_vec = [];
    a = 1;
    b = length(exp_vec);
    amp_diff = exp_vec(b)-(exp_vec(b) - exp_vec(a))/2;
    len_diff = distance_vec(b) - distance_vec(a);
    area_T = amp_diff*len_diff;

    amp_diff_A = exp_vec(fix(b/2)) - (exp_vec(fix(b/2)) - exp_vec(a))/2;
    len_diff_A = distance_vec(fix(b/2)) - distance_vec(a);
    area_A = amp_diff_A*len_diff_A;

    amp_diff_B = exp_vec(b) - (exp_vec(b) - exp_vec(fix(b/2)))/2;
    len_diff_B = distance_vec(b) - distance_vec(fix(b/2));
    area_B = amp_diff_B*len_diff_B;
    
    error = 0.1*abs(area_T-area_A-area_B);
    eval = error < tar;
    
    if error < tar
        fprintf('a ind = %i, b ind = %i, ', a, b)
        fprintf('a = %0.4i, b = %0.4i\n', distance_vec(a),distance_vec(b))
        end_vec = distance_vec(b);
    else
        right = simpsonFunc(exp_vec(fix(end/2):end),distance_vec(fix(end/2):end),tar);
        left = simpsonFunc(exp_vec(1:fix(end/2)),distance_vec(1:fix(end/2)),tar); 
        if isempty(end_vec)
            end_vec = [left right];
        else
            end_vec = [end_vec left right];
        end
    end
end