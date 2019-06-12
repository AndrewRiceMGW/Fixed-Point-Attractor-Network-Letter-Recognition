%% Simulator - Recurrent (Attractor) Network - Continuous Time
%% D. Bernal-Casas
%% Version 10: 2019 - Adapted from the literature
%% Main Sources:
%% [1] Trappenberg 
%% [2] Dayan and Abbot
%% [3] Churchland and Sejnowski
%% [4] Rolls and Deco
%% [5] Izhikevich
%% [6] Wilson
%% [7] Gerstner

clear all;
close all;
hold off; 
clc;

%% parameters of the model

nn=156; dx=1/nn; 
% We are setting C equals to 0 because we already take into account
% inhibition (negative values) in the pattern formation
C=0; 
dt =0.1;
%% Training weight matrix

%% The overlap here is
%% the normalized dot product of the network
%% states during an update with all
%% of the 10 patterns that were imprinted
%% with Hebbian learning into the network.

rng('default');
rng(5);
% pat = floor(2*rand(nn, 10))-0.5;
load pattern1;
pat = reshape(pattern1', 26, nn);
pat=floor(2*pat')-0.5;
w=pat*pat'; 
w=w/w(1,1);
w=100*(w-C);
%% Update with localised input 
tall = []; rall = [];

%% A noisy version of one stored pattern was
%% applied as external input until t = 10?.
 
%% Selects the first pattern
% This is the noisy version of the input pattern number 1   
%% A. I_ext=pat(:,1)+0.5;
% I_ext=pat(:,1)+0.5;
I_ext = zeros(1, nn)'; 
%% I_ext=pat(:,1);
% Number of neurons that we perturb - we can play with this number
% For a given number of neurons the input pattern is not correct
% We flip the responses
% I_ext(1:10)=1-I_ext(1:10);
[t,u]=ode45('raf',[0 20],zeros(1,nn),[],nn,dx,w,I_ext);
r=u>0.; tall=[tall;t]; rall=[rall;r];


%% Update without input 
% I_ext=zeros(nn,1);
Error = [];
for i = 1:15
    pat2 = pattern1 - (i*0.20)* rand(size(pattern1));
%     pat2 = imnoise(pattern1,'salt & pepper', i*0.00001);
    pat2 = reshape(pat2', 26, nn);
%     imshow(pat2);
    
    I_ext2 = pat2(1,:)';
%     I_ext2(1:10,:) = 1-I_ext(1:10,:);
% We set the initial condition to u(size(u,1),:) in order to have
% continuous differential equations
    [t,u]=ode45('raf',[0 20],zeros(1,nn),[],nn,dx,w,I_ext2);
    r_1=u>0.; tall=[tall;t]; rall=[rall;r_1];
    
    Error(i) = 0.5*sum(sum(r_1 - r))^2/nn;  
end
%% Plotting results
figure(2) 
plot(Error)
xlabel('Noise'); ylabel('Error'); 
%% Normalized dot product of the network
% %% A. plot(tall,4*(rall-0.5)*pat/nn)
% figure(3)
% plot(tall, rall)
% % plot(tall,4*(rall-0.5)*pat/nn)
% % % plot(tall,4*(rall)*pat/nn)
% % legend('1','2','3','4','5','6','7','8','9','10')