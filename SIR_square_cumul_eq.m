% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------




function dy = SIR_square_cumul_eq(t,y)

global beta alphai C alphai_min gamma epsi Dmaxi Dstar choice


S = y(1);
I = y(2);
R = y(3);
Cm = y(4);
dy = zeros(4,1);

if choice ==1
beta = C*(alphai_min + (alphai-alphai_min)/2*(1-tanh((Dstar-Dmaxi)*epsi)));
end
if choice == 2
%linear function alpha(D)
if Dstar<=Dmaxi*2
beta = C*(-(alphai-alphai_min)/(Dmaxi*2)*Dstar+alphai);
elseif Dstar>Dmaxi*2
beta = C*alphai_min;
end
end





dy(1) = -beta*I*S;
dy(2) = beta*I*S-gamma*I;
dy(3) = gamma*I;
dy(4) = beta*I*S;

end