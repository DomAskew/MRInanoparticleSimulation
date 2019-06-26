%
% 3D Walk - NANO-particle present
%
clear all; clc; clf;
%% SET CONSTANTS:
N = 10000; % Number of particles.
F = 2000; % Frequency of steps /Hz.
dt = 1 / F; % Computing time step /s.
t = 3; % Duration of walk /s.
n =( t / dt)+1; % Calculating total steps.
D = 3.08e-14; % Diffusion constant /m^2/s.
Bound = 1.87e-6; % Boundary condition.
% boundary condition to be 1Mol/m^3 being volume of L = 2.5e-7m.
%% RANDOM WALK 3D:
WalkX = zeros(N,n);
WalkX1 = zeros(N,n);
WalkX1(1:N,1) = linspace(Bound,-Bound,N);
WalkY = zeros(N,n);
WalkY1 = zeros(N,n);
WalkY1(1:N,1) = linspace(Bound,-Bound,N);
WalkZ = zeros(N,n);
WalkZ1 = zeros(N,n);
WalkZ1(1:N,1) = linspace(Bound,-Bound,N);
for I = 1:N % X - component walk.
    WalkX(I,:) = Diffusionwalk1(n, dt, D, Bound, WalkX1(I,:));
end
for J = 1:N % Y - component walk.
    WalkY(J,:) = Diffusionwalk2(n, dt, D, Bound, WalkY1(J,:));
end
for K = 1:N % Z - component walk.
    WalkZ(K,:) = Diffusionwalk3(n, dt, D, Bound, WalkZ1(K,:));
end
WalkR = sqrt( (WalkX.^2) + (WalkY.^2) + (WalkZ.^2) ); % Radial Component.
WalkT = acos(WalkZ ./ WalkR); % Angle to Z-axis /RAD.
%% NANOPARTICLE:
Mx = 6e6; % Magnetic Susceptibility.
H = 9.4; % External field strength /T.
ND = 30e-9; % Diameter of Nanoparticle /m.
MAG = Mx * H; % Magnetization /A/m.
Volume = 1/6 * pi * ND^3; % Volume of Nanoparticle /m^3.
% MM = MAG * Volume; % 'Induced' Magnetization is defined as the dipole moment per unit volume /Am^2.
MM = 5e-16;
u0 = 1.257e-6; % Permeability of free space /mkg/s^2A^2.
Bz = (( u0 * MM )/( 4 * pi )) * ((( 3 .* (cos(WalkT).^2) ) - 1 )./( WalkR.^3 )); % Magnetic field /T.
Gamma = 267.513e6; % Gyromagnetic ratio /RAD/sT.
%% NANOPARTICLE B FIELD:
[Y,Z] = meshgrid(linspace(-1e-5,1e-5,500));
R = sqrt( (Y.^2) + (Z.^2) );
Thetacos = Z ./ R;
Bz1 = (( u0 * MM )/( 4 * pi )) * ((( 3 .* (Thetacos.^2) ) - 1 )./( R.^3 ));
contour3(Y,Z,Bz1,230)
hold on
% axis([ -1.5e-7 1.5e-7 -1.8e-7 1.8e-7])

% Plot single walk
%stem3(WalkX(1,:),WalkY(1,:),WalkZ(1,:))

%% Spin Decay:
figure
Walklarmor = Gamma .* Bz;
Walkphase = (dt .* Walklarmor(:,1:n)); % relative phase at each point due to larmor frequency.
for L = 1:N % Loop to calculate progressive phase per particle.
    Walkphase(L,:) = cumsum(Walkphase(L,:));
end
XPHASE = cos(Walkphase); % X component of phase.
YPHASE = sin(Walkphase); % Y component of phase.
AVGX = mean(XPHASE); 
AVGY = mean(YPHASE);
MAGXY = ( (AVGX.^2) + (AVGY.^2) ).^0.5; % Pythagoras to determine Magnetization.
timestep = [1:n] .* dt; % Time step array for plotting.
plot(timestep,MAGXY) % plot of Magnetization decay.

%% Spin Echo
figure
Ninitial = 0.01 / dt; % initial positioning of particles.
noff = 0.15 / dt; % step number at which field is turned off.
non = 0.15 / dt; % step number at which field is turned back on after pi pulse.
npulse = ( non + noff )/2;
noff1 = 1.50 / dt; % step number at which field is turned off again.
Walklarmor1 = zeros(N,n); %
XPHASE1 = zeros(N,n); % setup array of dimentions N,n.
YPHASE1 = zeros(N,n); %
Walklarmor1(:,Ninitial:noff) = Gamma .* Bz(:,Ninitial:noff); % larmor frequency like before, during the period of active B field.
% between noff and non, the field gradient g = 0 therefore larmor frequency
% is equal to 0.    
Walklarmor1(:,(non+1):noff1) =  Gamma .* Bz(:,(non+1):noff1); % larmor frequency when G=1 again after pi pulse.
%after noff1, the field gradient G = 0 again therefore larmor freq = 0.
Walkphase1 = (dt .* Walklarmor1); % relative phase at each point.
for L1 = 1:N % Loop to find the progressive phase change per particle over time up to the pi pulse.
    Walkphase1(L1,1:npulse) = cumsum(Walkphase1(L1,1:npulse));
end
XPHASE1(:,1:npulse) = cos(Walkphase1(:,1:npulse));
YPHASE1(:,1:npulse) = sin(Walkphase1(:,1:npulse));
Walkphase1(:,npulse) = -1 * Walkphase1(:,npulse); % Pi pulse - change sign.
for L2 = 1:N % Loop to find the progressive phase change per particle over time after pi pulse.
    Walkphase1(L2,npulse:n) = cumsum(Walkphase1(L2,npulse:n));
end
XPHASE1(:,npulse:n) = cos(Walkphase1(:,npulse:n));
YPHASE1(:,npulse:n) = sin(Walkphase1(:,npulse:n));
AVGX1 = mean(XPHASE1);
AVGY1 = mean(YPHASE1);
MAGXY1 = ( (AVGX1.^2) + (AVGY1.^2) ).^0.5;
timestep = [1:n] .* dt; % Time step array for plotting.

plot(timestep,MAGXY1) % plot of this new average phase against time.
xlabel('Time (s)')
ylabel('Magnetization (A/m)')
axis([0 1 0 1])
%%
% magnetic moment variation:
Mmoment = [1e-16, 3e-16, 6e-16, 8e-16, 1e-15];%values read from simulations.
trelax =1e3.* 10./[38.99, 63.69, 96.01, 107.52, 117.98];
trelax1 =1e3.* 10./[38.99, 63.69, 86.01, 105.52, 120.98];
plot(Mmoment, trelax1,'r*')
xlabel('Magnetic moment (Am^2)')
ylabel('T_2 Relaxation (ms)')
hold on
fitobject=linspace(min(Mmoment), max(Mmoment));
fit1 = interp1(Mmoment, trelax, fitobject, 'spline');
plot(fitobject, fit1)
 
%%
% Concentration variation:
ConcentrationT2 = [1, 2.5, 4, 6, 7.5]; %values read from simulations.
trelax2 = 1e3.*[0.238, 0.127, 0.085, 0.065, 0.057];
plot(ConcentrationT2,trelax2,'r*')
xlabel('Concentration (mmolm^-^3)')
ylabel('T_2 Relaxation (ms)')
hold on
fitobject2=linspace(min(ConcentrationT2), max(ConcentrationT2));
fit12 = interp1(ConcentrationT2, trelax2, fitobject2, 'spline');
plot(fitobject2, fit12)
axis([1 7.5 40 240])