% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------




function dy = onlyD_eq_updated(t,y)

global C delta Hmax k h


H = y(1);
dy = zeros(1,1);


dy(1) = C*k*h*1/2*(1-tanh(3*(H-Hmax))) - delta*H;

end