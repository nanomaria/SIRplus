% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------


clear all
close all


alpha_max = 0.04;
alpha_min = 0.02;

p_max = 0.05;
p_min = 0.03;

Hmaxi = 5; 
Hmaxm = 10; 


D = 0:.1:15;



alpha_lin = zeros(1,length(D));
for i = 1:size(D,2)
if D(i) <= Hmaxi*2
alpha_lin(i) = -(alpha_max-alpha_min)/(Hmaxi*2)*D(i)+alpha_max;
elseif D(i)>Hmaxi*2
alpha_lin(i) = alpha_min;
end
end


figure(3)
subplot(1,2,1)
plot(D,(alpha_min + (alpha_max-alpha_min)/2*(1-tanh((D-Hmaxi)*1))), 'Color',[0.995 0.681 0],'Linestyle','-', 'Linewidth',1.8);
hold on
plot(D,alpha_lin,'Color',[0.995 0.681 0],'Linestyle','--', 'Linewidth',1.8);
axis([0  11 0.015 alpha_max*(1.1)])
legend({'Sigmoid','Linear'},'FontSize',14)
xlabel('Interaction-dependent health benefits (H)')
ylabel('\alpha(H)')
set(gca,'fontsize', 14)


p_lin = zeros(1,length(D));
for i = 1:size(D,2)
if D(i) <= Hmaxi*2
p_lin(i) = -(p_max-p_min)/(Hmaxi*2)*D(i)+p_max;
elseif D(i)>Hmaxi*2
p_lin(i) = p_min;
end
end



figure(3)
subplot(1,2,2)
plot(D,(p_min + (p_max-p_min)/2*(1-tanh((D-Hmaxi)*1))), 'Color',[0.592 0.054 0.325],'Linestyle','-', 'Linewidth',1.8);
hold on
plot(D,p_lin,'Color',[0.592 0.054 0.325],'Linestyle','--', 'Linewidth',1.8);
axis([0  11 0.025 p_max*(1.1)])
legend({'Sigmoid','Linear'},'FontSize',14)
xlabel('Interaction-dependent health benefits (H)')
ylabel('p(H)')
set(gca,'fontsize', 14)

