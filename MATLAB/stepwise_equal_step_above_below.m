%Stepwise approximation of exponential DDF

%%
clc; clear; close all;

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

tot_steps = length(z_vector); %why 63?
segs = 4; %input variable
b2_step = tot_steps/segs;
b2_vec = b2_initial;
b2_vec_cont = b2_initial*exp(-sigma*z_vector);
b2_vec2 = b2_initial*exp(-0.08*z_vector);
b2_vec3 = b2_initial*exp(-0.01*z_vector);
b2_vec4 = b2_initial*exp(-0.1*z_vector);
%%
figure
semilogy(z_vector, abs(b2_vec_cont))
hold on
semilogy(z_vector,abs(b2_vec2))
hold on
semilogy(z_vector,abs(b2_vec3))
hold on
semilogy(z_vector,abs(b2_vec4))
legend('0.05','0.08','0.01', '0.09');
ylim([7.4 33])
%%
i_vec = 1:length(z_vector);
i_sel_vec = i_vec(fix(mod(i_vec,tot_steps/segs)) == 0);
for i = 1:3
    i_half_vec(i) = fix(i_sel_vec(i) + (i_sel_vec(i)-i_sel_vec(i+1))/2);
end
%i_half_vec(4) = fix(i_sel_vec(4) - (i_sel_vec(4)-i_sel_vec(3))/2);


% i_half_vec = [123 248 372 497];

for i = 1:4
    b2_half_vec(i) = b2_vec_cont(i_sel_vec(i));
end


b2_track = 1;
%%
for i = 1:length(z_vector)


    D = exp(1i*(dz/2)*b2*f.^2).*fft(A);
    D = ifft(D);
    N = exp(1i*gama*abs(A).^2*dz);
    A = D.*N;
    spec_wave(i,:) = ifftshift(abs(fft(A)));
    simul_wave(i,:) = abs(A).^2;
    

%     plot(t/To, simul_wave(i,:));
%     ylim([0 11])
%     pause(.00000001);
        
%         if fix(mod((i),tot_steps/segs)) == 0
        if find(i_half_vec == i) > 0
%             b2 = b2_initial*exp(-sigma*z_vector(i)); 
            b2 = b2_half_vec(b2_track);
            b2_track = b2_track + 1;
            b2_vec = [b2_vec b2];
        else
            b2_vec = [b2_vec b2];
        end
        
    
        
    fwhm=find(abs(abs(A).^2)>abs(max(abs(A).^2)/2));
    fwhm=length(fwhm);
    fwhm_vec(i) = fwhm;
    c_factor(i) = (fwhm1/fwhm);
    
%     if i == 885
%         break;
%     end
    

end

%%
figure;
mesh(t/To,z_vector, simul_wave);
xlabel('Time Delay T/To');
ylabel('Distance z/LD');
zlabel('Amplitude');
title('First order soliton with initial chirp');

%%
figure;
plot(c_factor,'-o','MarkerIndices',i_half_vec,'MarkerFaceColor','red')
title('compression factor')
max_cf = max(c_factor);
%%
% target = 987;
% eT = trapz(t/To,simul_wave(target,:));
% simAmp = max(simul_wave(target,:));
% finWidth = fwhm_vec(target);
% eSech = 2*simAmp*((finWidth*dt/To)/1.763); %1.763 is a constant from Cao Wai 2005 paper
% PE = abs(eT-eSech)/eT *100;
%%
for i = 1:length(z_vector)
    eT(i) = trapz(t/To,simul_wave(i,:));
    simAmp(i) = max(simul_wave(i,:));
    eSech(i) = 2*simAmp(i)*((fwhm_vec(i)*dt/To)/1.763); %1.763 is a constant from Cao Wai 2005 paper
    PE(i) = abs(eT(i)-eSech(i))/eT(i) *100;
end

% finWidth = fwhm_vec(target);
figure
plot(PE)
hold on
plot(c_factor)
% plot(normalize(PE,'range'))
% hold on;
% plot(normalize(c_factor,'range'))
% hold on;
% plot(normalize(c_factor,'range')-normalize(PE,'range'))


%%
b2_vec = [b2_vec ones(1,length(z_vector)-length(b2_vec))*b2];
%%
figure
plot(z_vector,b2_vec(1:end-1))
hold on;
plot(z_vector, b2_vec_cont)
% title('above and below')
xlabel('Distance z');
ylabel('Dispersion coefficient \beta_2');
% xlim([-30 0])
% hold on;
% i_vec_2 = [i_vec(1:886), fliplr(i_vec(1:886))];
% inBetween = [b2_vec(1:886), fliplr(b2_vec_cont(1:886))];
% fill(i_vec_2, inBetween, 'g');
%%
area_diff = trapz(abs(b2_vec)) - trapz(abs(b2_vec_cont));
%%
% figure
% plot(simul_wave(target,:));
% hold on;
% plot(simul_wave(target-1,:))
% hold on;
% plot(simul_wave(target+1,:))
% legend('tar','tar-1','tar+1');