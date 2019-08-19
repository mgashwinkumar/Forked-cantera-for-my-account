function [T,rho,x,Y]=PFR_setup(A_in,A_out,L,n,k,gas,mdot)

dAdx=abs(A_in-A_out)/L;
dx=L/n;% The whole length of the reactor is divided into n small lengths

T(1)=temperature(gas);
P(1)=pressure(gas);
Y(1,:)=massFractions(gas);    
rho(1)=density(gas);

x=0:dx:L;

for i=2:length(x)
    
    i
    Total_i=length(x)
    nsp = nSpecies(gas);
    initial(1) = rho(i-1);       %---------------------------------------------------------------------------------------------------------
    initial(2) = T(i-1);         %----The values of variables at previous location are given as initial values to the current iteration----
    initial(3:nsp+2) = Y(i-1,:); %---------------------------------------------------------------------------------------------------------
    limits=[x(i-1),x(i)];        % The limmits of the current reactor
    set(gas,'T',T(i-1),'Density',rho(i-1),'MoleFractions',Y(i-1,:));
    options = odeset('RelTol',1.e-10,'AbsTol',1e-10,'InitialStep',1e-8,'NonNegative',1);
    [x_soln,y]= ode15s(@PFR_solver_CANTERA,limits,initial,options,gas,mdot,A_in,dAdx,k); % These values are passed onto the ode15s solver
   
    T(i)=y(end,2);          %----------------------------------------------
    rho(i)=y(end,1);        %-----The last location values are copied------
    Y(i,:)=y(end,3:nsp+2);  %----------------------------------------------

end



end
