% This code snippet is a sample to model a converging nozzle as
% Plug Flow Reactor with given dimensions and an incoming gas. The pressure
% is not constant here as opposed to the Python Version of the Plug Flow
% Reactor modelling.

%% Setting the Gas
clear all
close all
clc
gas_conv = IdealGasMix('gri30.cti');
ich4 = speciesIndex(gas_conv,'CH4');
io2  = speciesIndex(gas_conv,'O2');
in2  = speciesIndex(gas_conv,'N2');
nsp = nSpecies(gas_conv);
x = zeros(nsp,1);
% Change the below values for different Phi values of methane Combsution
x(ich4,1) =0.2899; 
x(io2,1) = 2.0;
x(in2,1) = 7.52;
set(gas_conv,'T',1473,'P',4.4763e5,'MoleFractions',x);
gas_conv=equilibrate(gas_conv,'HP');

%% Converging Section
% The Dimensions and conditions of the reactor are given below
A_in=0.018;            % Inlet Area
A_out=0.003;           % Exit Area
L=1.284*0.0254;        % Length of the reactor
n=100;                 % The whole reactor is divided into n small reactors
mdot_conv=1.1466;     % Mass flow rate into the reactor
% k = -1 makes the solver solve for converging area. 
% Use k = +1 for diverging area. make sure to change the inlet and exit area for diverging area
% Use k = 0 for constant cross sectional area
k=-1;     
%--------------------------------------------------------------
dAdx=abs(A_in-A_out)/L;
[T_conv,rho_conv,x_conv,Y_conv]=PFR_setup(A_in,A_out,L,n,k,gas_conv,mdot_conv); % The above values are passed on to the PFR function to solve for density, temperature and Mas fractions

for i=1:length(x_conv)
    gas_after_solving_conv = gas_conv;
    set(gas_after_solving_conv,'Temperature',T_conv(i),'Density',rho_conv(i),'MassFractions',Y_conv(i,:)); % The gas is set to the solved property values at each location
    A_conv(i) = A_in+k*x_conv(i)*dAdx;
    vx_conv(i) = mdot_conv./(A_conv(i)*rho_conv(i)); % Velocity is calculated from Mass flow rate, Area and Density
    R_conv(i)=8314/meanMolecularWeight(gas_after_solving_conv); % Specific Gas Constant
    M_conv(i) = vx_conv(i)/soundspeed(gas_after_solving_conv); % Mach  No. is calculated from local velocity and local speed of sound
    P_conv(i) = rho_conv(i)*R_conv(i)*T_conv(i); % Pressure is calculated from density, temeprature adn gas constant
    disp('Finished Solving Convergent part, setting gas now')
    i
    length(x_conv)
end

%% Plotting
plot(x_conv,M_conv)
xlabel('Position')
ylabel('Mach No')
title('Mach No Variation')
figure(2)
plot(x_conv,A_conv)
xlabel('Position')
ylabel('Area(m^2)')
title('Reactor Profile')
figure(3)
plot(x_conv,T_conv)
xlabel('Position')
ylabel('Temperature')
title('Temperature Variation')
figure(4)
plot(x_conv,rho_conv)
xlabel('X-Position (m)')
ylabel('Density (Kg/m^3)')
title('Density Variation')
figure(5)
plot(x_conv,P_conv)
xlabel('X-Position (m)')
ylabel('Pressure (Pa)')
title('Pressure Variation')
