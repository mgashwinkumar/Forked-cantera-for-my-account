function [F] = PFR_solver_CANTERA(x,initial,gas,mdot,A_in,dAdx,k)

rho = initial(1);
T = initial(2);
Y = initial(3:length(initial));

if k==1
    A = A_in-k*dAdx*x;
elseif k==-1|k==0
    A=A_in+k*dAdx*x;
    dAdx=-dAdx;
end

MW_mix = meanMolecularWeight(gas);
Ru = 8314;
R = Ru/MW_mix;
nsp = nSpecies(gas);
vx = mdot/(rho*A);
P = rho*R*T;
if k~=0
    set(gas,'Temperature',T,'Density',rho,'MassFractions',Y); % the gas is set to the corresponding properties during each iteration of the ode loop
else
    set(gas,'Temperature',T,'Pressure',P,'MassFractions',Y); % the gas is set to the corresponding properties during each iteration of the ode loop

end
MW = molecularWeights(gas);
h = enthalpies_RT(gas).*R.*T;
w = netProdRates(gas);
%--------------------------------------------------------------------------
%---F(1), F(2) and F(3:end) are the differential equations modelling the---
%---density, temperature and mass fractions variations along a plug flow---
%-------------------------reactor------------------------------------------
%--------------------------------------------------------------------------
F(1)=((1-R/cp_mass(gas))*((rho*vx)^2)*(1/A)*(dAdx) + rho*R*sum(MW.*w.*(h-MW_mix*cp_mass(gas)*T./MW))/(vx*cp_mass(gas)) ) / (P*(1+vx^2/(cp_mass(gas)*T)) - rho*vx^2);
F(2) = (vx*vx/(rho*cp_mass(gas)))*F(1) + vx*vx*(1/A)*(dAdx)/cp_mass(gas) - (1/(vx*rho*cp_mass(gas)))*sum(h.*w.*MW);

for i = 3:nsp+2
    F(i) = w(i-2)*MW(i-2)/(rho*vx);
end
F = F';

end

