% SIR+ models: Accounting for interaction-dependent disease susceptibility in the planning of public health interventions
%------------------------------------------------------------------------
% AUTHORSHIP AND CONTACT
% Original code by Maria M. Martignoni (Hebrew University of Jerusalem)
% Contact: maria.martignonimseya@mail.huji.ac.il
%-------------------------------------------------




function dy = onlyD_eq(t,y)

global C delta epsm Dmaxm k sR iR rR pi pr 


D = y(1);
dy = zeros(1,1);

        
%zeta = C*k*1/2*(1-tanh(epsm*(D-Dmaxm)));

%the 0.5 is for the tanh function!
dy(1) = C*k*1/2*(1-tanh(epsm*(D-Dmaxm)))*(sR+pi*iR+pr*rR)^2 - delta*D;

end