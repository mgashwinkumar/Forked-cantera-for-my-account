% General_Plug_Flow_Reactor(PFR) - to solve PFR equations for reactors
%
%    This code snippet is to model a constant area and varying area
%    (converging and diverging) nozzle as Plug Flow Reactor with given 
%    dimensions and an incoming gas. The pressure is not assumed to be
%    constant here as opposed to the Python Version of the
%    Plug Flow Reactor modelling.
%    
%    The governing equations used in this code can be referenced at:
%    *S.R Turns, An Introduction to Combustion - Concepts and Applications,
%    McGraw Hill Education, India, 2012, 206-210.*
%

%% Setting the Gas
clear all
close all
clc

%Temperature of gas
T0=1473;
%Pressure of gas
P0=4.47*101325;

%Equivalence Ratio
Phi=0.2899;

gas_calc = IdealGasMix('gri30.cti');
ich4 = speciesIndex(gas_calc,'CH4');
io2  = speciesIndex(gas_calc,'O2');
in2  = speciesIndex(gas_calc,'N2');
nsp = nSpecies(gas_calc);
x = zeros(nsp,1);
% Change the below values for different Phi values of methane Combustion
x(ich4,1) = Phi; 
x(io2,1) = 2.0;
x(in2,1) = 7.52;
set(gas_calc,'T',T0,'P',P0,'MoleFractions',x);
gas_calc=equilibrate(gas_calc,'HP');

%% Calculation of properties along the reactor length
% The Dimensions and conditions of the reactor are given below

% Inlet Area
A_in = 0.018;
% Exit Area
A_out = 0.003;
% Length of the reactor
L = 1.284*0.0254;  
% The whole reactor is divided into n small reactors
n = 100;
% Mass flow rate into the reactor in Kg/s
mdot_calc = 1.125;
% k = -1 makes the solver solve for converging area. 
% k = +1 makes the solver solve for diverging area.
% k = 0 makes the solver solve for constant cross sectional area
if A_in>A_out
    k=-1;
elseif A_out>A_in
    k=1;
else k=0;
end
    
% The above values are passed on to the PFR function to solve for density, temperature and Mass fractions
dAdx = abs(A_in-A_out)/L;
[T_calc,rho_calc,x_calc,Y_calc] = PFR_setup(A_in,A_out,L,n,k,gas_calc,mdot_calc); 
A_calc = A_in+k.*x_calc*dAdx;
for i=1:length(x_calc)
    gas_after_solving_calc = gas_calc;
    % The gas is set to the solved property values at each location
    set(gas_after_solving_calc,'Temperature',T_calc(i),'Density',rho_calc(i),'MassFractions',Y_calc(i,:));
    % Velocity is calculated from Mass flow rate, Area and Density
    vx_calc(i) = mdot_calc./(A_calc(i)*rho_calc(i));
    % Specific Gas Constant
    R_calc(i) = 8314/meanMolecularWeight(gas_after_solving_calc);
    % Mach  No. is calculated from local velocity and local speed of sound
    M_calc(i) = vx_calc(i)/soundspeed(gas_after_solving_calc);
    % Pressure is calculated from density, temeprature and gas constant
    P_calc(i) = rho_calc(i)*R_calc(i)*T_calc(i);
end

%% Plotting
plot(x_calc,M_calc)
xlabel('X-Position (m)')
ylabel('Mach No')
title('Mach No Variation')
figure(2)
plot(x_calc,A_calc)
xlabel('X-Position (m)')
ylabel('Area(m^2)')
title('Reactor Profile')
figure(3)
plot(x_calc,T_calc)
xlabel('X-Position (m)')
ylabel('Temperature')
title('Temperature Variation')
figure(4)
plot(x_calc,rho_calc)
xlabel('X-Position (m)')
ylabel('Density (Kg/m^3)')
title('Density Variation')
figure(5)
plot(x_calc,P_calc)
xlabel('X-Position (m)')
ylabel('Pressure (Pa)')
title('Pressure Variation')
