clc, close all, clear all;

%This is the exponential function setup
b2_initial = -20;

To = 10;
C = -0.01;
ao = -4/((To^2) * pi);
% sigma = C*b2_initial;
sigma = ao*b2_initial;

To= 10; %Pulse width
Ld = To^2/abs(b2_initial);
Lnl = Ld;
dz = min(Ld,Lnl)/100;
z_vector = 0:dz:12;

b2_vec_cont = abs(b2_initial*exp(-sigma*z_vector));

% %%
% a = 1;
% b = length(b2_vec_cont);
% amp_diff = abs(b2_vec_cont(b) - b2_vec_cont(a));
% len_diff = z_vector(b) - z_vector(a);
% area = amp_diff*len_diff;
% 
% amp_diff_A = abs(b2_vec_cont(b/2) - b2_vec_cont(a));
% len_diff_A = z_vector(b/2) - z_vector(a);
% area_A = amp_diff_A*len_diff_A;
% 
% amp_diff_B = abs(b2_vec_cont(b) - b2_vec_cont(b/2));
% len_diff_B = z_vector(b) - z_vector(b/2);
% area_B = amp_diff_B*len_diff_B;
% 
% error = 0.1*abs(area-area_A-area_B);

%%
vec = simpsonFunc(b2_vec_cont, z_vector,0.15);
%%
for i = 1:length(vec)
    vec(i) = find(z_vector == vec(i));
end
vec = [0 vec length(z_vector)];
% %%
% seg_stops = [0];
% 
% % seg_stops = [seg_stops find(z_vector == 3.6667e+00)];
% % seg_stops = [seg_stops find(z_vector == 7.4242e+00)];
% % seg_stops = [seg_stops find(z_vector == 1.4939e+01)];
% % seg_stops = [seg_stops find(z_vector == 2.2424e+01)];
% % seg_stops = [seg_stops find(z_vector == 0030)];
% 
% % seg_stops = [seg_stops 246 494 741 991];
% seg_stops = [60 122 246 494 991];
% vec = [0 vec 991];
%%
% seg_vec = [];
% for i = 1:length(seg_stops)-1
%    segment_len =  seg_stops(i+1) - seg_stops(i);
%    segment = ones(1,fix(segment_len))*b2_vec_cont(fix(segment_len/2 + seg_stops(i)));
%    seg_vec = [seg_vec segment];
% end
% seg_vec = [seg_vec ones(1,length(z_vector)-length(seg_vec))*seg_vec(end)];

seg_vec2 = [];
for i = 1:length(vec)-1
   segment_len =  vec(i+1) - vec(i);
   segment = ones(1,fix(segment_len))*b2_vec_cont(fix(segment_len/2 + vec(i)));
   seg_vec2 = [seg_vec2 segment];
end
seg_vec2 = [seg_vec2 ones(1,length(z_vector)-length(seg_vec2))*seg_vec2(end)];
% %%
% figure
% plot(z_vector,(-1)*seg_vec);
% hold on
% plot(z_vector,(-1)*b2_vec_cont);
% % title('Adaptive Midpoint rule');
% xlabel('Distance z');
% ylabel('Dispersion coefficient \beta_2');

%%
figure
plot(z_vector,(-1)*seg_vec2);
hold on
plot(z_vector,(-1)*b2_vec_cont);
title('Automated');
xlabel('Distance z');
ylabel('Dispersion coefficient \beta_2');
%%
% function area = simpsonFunc(exp_vec, distance_vec,tar)
% 
% %     figure
% %     plot(distance_vec,exp_vec)
% %     trapz(distance_vec,exp_vec);
%     a = 1;
%     b = length(exp_vec);
%     amp_diff = exp_vec(b)-(exp_vec(b) - exp_vec(a))/2;
%     len_diff = distance_vec(b) - distance_vec(a);
%     area_T = amp_diff*len_diff;
% 
%     amp_diff_A = exp_vec(fix(b/2)) - (exp_vec(fix(b/2)) - exp_vec(a))/2;
%     len_diff_A = distance_vec(fix(b/2)) - distance_vec(a);
%     area_A = amp_diff_A*len_diff_A;
% 
%     amp_diff_B = exp_vec(b) - (exp_vec(b) - exp_vec(fix(b/2)))/2;
%     len_diff_B = distance_vec(b) - distance_vec(fix(b/2));
%     area_B = amp_diff_B*len_diff_B;
%     
%     error = 0.1*abs(area_T-area_A-area_B);
%     eval = error < tar;
%     
%     if error < tar
%         fprintf('tar %0.2i, error %0.2i, eval %i, ', tar,error,eval)
%         fprintf('a = %0.4i, b = %0.4i\n', distance_vec(a),distance_vec(b))
%         area = abs(area_T);
%     else
%         right = simpsonFunc(exp_vec(fix(end/2):end),distance_vec(fix(end/2):end),tar);
%         left = simpsonFunc(exp_vec(1:fix(end/2-1)),distance_vec(1:fix(end/2-1)),tar); 
%         area = left+right;
%     end
% end

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