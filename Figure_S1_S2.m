% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------



close all

C= 5;

alphai_max = 0.1;
alphai_min = 0.02;

pd_max = 0.1;
pd_min = 0.02;

epsi = 1;   
epsm = 1;
Dmaxi = 3; %if we have this IDHB diversity, it's hard to get more
Dmaxm = 6; % after this, more diversity doesn't help to prevent more infection


D = 0:.1:15;

T = zeros(1,length(D));
for i = 1:size(D,2)
if D(i) <= Dmaxm*2
T(i) = -1/(Dmaxm*2)*D(i)+1;
elseif D(i)>Dmaxm*2
T(i) = 0;
end
end

figure(4)
plot(D,0.5*(1-tanh(10*(D-Dmaxm))), 'k-', 'Linewidth',2.2);
hold on
plot(D,0.5*(1-tanh(1*(D-Dmaxm))), 'k-', 'Linewidth',1.8);
hold on
plot(D,0.5*(1-tanh(0.5*(D-Dmaxm))), 'k-', 'Linewidth',1.4);
hold on
plot(D,T,'k-', 'Linewidth',1);
legend('\epsilon_m = 10','\epsilon_m = 1','\epsilon_m = 0.5', 'Linear')
axis([0  D(end) 0 (1+0.1)])

xlabel('Interaction-dependent health benefits (H)')
ylabel('T(H)')
set(gca,'fontsize', 14)



alpha_lin = zeros(1,length(D));
%linear function T(D)
for i = 1:size(D,2)
if D(i) <= Dmaxi*2
alpha_lin(i) = -(alphai_max-alphai_min)/(Dmaxi*2)*D(i)+alphai_max;
elseif D(i)>Dmaxi*2
alpha_lin(i) = alphai_min;
end
end


figure(3)
subplot(1,2,1)
plot(D,(alphai_min + (alphai_max-alphai_min)/2*(1-tanh((D-Dmaxi)*10))), 'Color',[0.995 0.681 0],'Linestyle','-', 'Linewidth',2.2);
hold on
plot(D,(alphai_min + (alphai_max-alphai_min)/2*(1-tanh((D-Dmaxi)*1))), 'Color',[0.995 0.681 0],'Linestyle','-', 'Linewidth',1.8);
hold on
plot(D,(alphai_min + (alphai_max-alphai_min)/2*(1-tanh((D-Dmaxi)*0.75))), 'Color',[0.995 0.681 0],'Linestyle','-', 'Linewidth',1.4);
hold on
plot(D,alpha_lin,'Color',[0.995 0.681 0],'Linestyle','-', 'Linewidth',1);
axis([0  8 0 alphai_max*(1.1)])
legend('\epsilon_i = 10','\epsilon_i=1','\epsilon_i = 0.75', 'Linear')
xlabel('Interaction-dependent health benefits (H)')
ylabel('\alpha(H)')
set(gca,'fontsize', 14)


pd_lin = zeros(1,length(D));
for i = 1:size(D,2)
if D(i) <= Dmaxi*2
pd_lin(i) = -(pd_max-pd_min)/(Dmaxi*2)*D(i)+pd_max;
elseif D(i)>Dmaxi*2
pd_lin(i) = pd_min;
end
end



figure(3)
subplot(1,2,2)
plot(D,(pd_min + (pd_max-pd_min)/2*(1-tanh((D-Dmaxi)*10))), 'Color',[0.592 0.054 0.325],'Linestyle','-', 'Linewidth',2.2);
hold on
plot(D,(pd_min + (pd_max-pd_min)/2*(1-tanh((D-Dmaxi)*1))), 'Color',[0.592 0.054 0.325],'Linestyle','-', 'Linewidth',1.8);
hold on
plot(D,(pd_min + (pd_max-pd_min)/2*(1-tanh((D-Dmaxi)*0.75))), 'Color',[0.592 0.054 0.325],'Linestyle','-', 'Linewidth',1.4);
hold on
plot(D,pd_lin,'Color',[0.592 0.054 0.325],'Linestyle','-', 'Linewidth',1);
axis([0  8 0 pd_max*(1.1)])
legend('\epsilon_i = 10','\epsilon_i=1','\epsilon_i = 0.75', 'Linear')
xlabel('Interaction-dependent health benefits (H)')
ylabel('p(H)')
set(gca,'fontsize', 14)
