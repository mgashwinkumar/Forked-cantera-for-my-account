function [T,rho,x,Y] = PFR_setup(A_in,A_out,L,n,k,gas,mdot)

dAdx = abs(A_in-A_out)/L;
dx = L/n;% The whole length of the reactor is divided into n small lengths

T(1) = temperature(gas);
P(1) = pressure(gas);
Y(1,:) = massFractions(gas);    
rho(1) = density(gas);

x = 0:dx:L;

for i = 2:length(x)
    
    %Solver location indicator
    Solving_current_reactor = i
    Total_reactors = length(x)
    nsp = nSpecies(gas);
%--------------------------------------------------------------------------
%------The values of variables at previous location are given as initial--- 
%------values to the current iteration and the limits of the current-------
%--------------reactor and the gas entering it are being set---------------
%--------------------------------------------------------------------------
    init(1) = rho(i-1);
    init(2) = T(i-1);
    init(3:nsp+2) = Y(i-1,:);
    limits = [x(i-1),x(i)];
    set(gas,'T',T(i-1),'Density',rho(i-1),'MoleFractions',Y(i-1,:));
    options = odeset('RelTol',1.e-10,'AbsTol',1e-10,'InitialStep',1e-8,'NonNegative',1);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
    % These values are passed onto the ode15s solver
    [~,y] = ode15s(@PFR_solver,limits,init,options,gas,mdot,A_in,dAdx,k);
   
    T(i) =y(end,2);
    rho(i)=y(end,1);
    Y(i,:)=y(end,3:nsp+2);
end



end
