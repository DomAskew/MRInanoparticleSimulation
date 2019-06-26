clear all; clc;
N = input('Number of particles: ');% set 10000.
dT = input('Set time step frequency: ');% steps per second.
dt = 1/dT; % computing timestep.
t = input('Duration of walk: ');% How long you wish the walk to run for/s.
n =( t / dt)+1; % finding total steps for arrays.
d = input('Diffusion constant: '); % 0.282cm^2/s.
bound = input('Set boundary condition(+-): ');
m=1; % 1D
Bzero = 9.4; % Default external B field.
G = 1;% Field gradient G/cm.
gamma = 42.57;%(2pi)(MHzT-1)
%
% Random walk for set particles N, from diffusion function.
%
walk = zeros(N,n);
dum = 1;
for dum = 1:N
    walk(dum,:) = Diffusion1Dwalk(m,dt,d,t,bound);
end
histpoints = walk(:,n);
histogram(histpoints) %plotting histogram of walk.
hold on
% Proof using formula:
x = linspace(-5,5,1001);
c = ((2*N*((pi*d)^0.5))/(2*9.1*((pi*d*t)^0.5))).*(exp((-1.*(x.^2))/(4*d*t)));
plot(x,c)
xlabel('Displacement (mm)')
ylabel('Number of Molecules')
%
% end
%
%
%  Larmor frequency at each point in random walk.
%
figure
larmorfrequency = (-1.*gamma).*Bzero; % larmor frequency.
G = 1.5; % Field gradient set to 1 g/cm.
walklarmor =  gamma .* walk .* G; % Larmor Frequency @ each point - constant B field.
walkphase = (dt .* walklarmor); % relative phase at each point due to larmor frequency.
for L = 1:N % Loop to calculate progressive phase per particle.
    walkphase(L,:) = cumsum(walkphase(L,:));
end
XPHASE = cos(walkphase); % X component of phase.
YPHASE = sin(walkphase); % Y component of phase.
AVGX = mean(XPHASE); 
AVGY = mean(YPHASE);
MAGXY = ( (AVGX.^2) + (AVGY.^2) ).^0.5; % Pythagoras to determine Magnetization.
timestep = [1:n] .* dt; % Time step array for plotting.
plot(timestep,MAGXY) % plot of Magnetization decay.
hold on
% Proof from Eq 3.50 of decay:
PHASEavg = exp((-1/3).*(gamma^2).*(G^2).*d.*(timestep.^3));
plot(timestep,PHASEavg)
xlabel('Time (s)')
ylabel('Magnetization (A/m)')
axis([0 0.7 0 1])
%
% end
%
%
% Spin echo
%
figure
Ninitial = 1.0 / dt; % initial positioning of particles.
noff = 1.07 / dt; % step number at which field is turned off.
non = 1.47 / dt; % step number at which field is turned back on after pi pulse.
npulse = ( non + noff )/2;
noff1 = 1.8 / dt; % step number at which field is turned off again.
walklarmor1 = zeros(N,n); %
XPHASE1 = zeros(N,n); % setup array of dimentions N,n.
YPHASE1 = zeros(N,n); %
walklarmor1(:,Ninitial:noff) = gamma .* walk(:,Ninitial:noff) .* G; % larmor frequency like before, during the period of active B field.
% between noff and non, the field gradient g = 0 therefore larmor frequency
% is equal to 0.    
walklarmor1(:,(non+1):noff1) =  gamma .* walk(:,(non+1):noff1) .* G; % larmor frequency when G=1 again after pi pulse.
%after noff1, the field gradient G = 0 again therefore larmor freq = 0.
walkphase1 = (dt .* walklarmor1); % relative phase at each point.
for L1 = 1:N % Loop to find the progressive phase change per particle over time up to the pi pulse.
    walkphase1(L1,1:npulse) = cumsum(walkphase1(L1,1:npulse));
end
XPHASE1(:,1:npulse) = cos(walkphase1(:,1:npulse));
YPHASE1(:,1:npulse) = sin(walkphase1(:,1:npulse));
walkphase1(:,npulse) = -1 * walkphase1(:,npulse); % Pi pulse - change sign.
for L2 = 1:N % Loop to find the progressive phase change per particle over time after pi pulse.
    walkphase1(L2,npulse:n) = cumsum(walkphase1(L2,npulse:n));
end
XPHASE1(:,npulse:n) = cos(walkphase1(:,npulse:n));
YPHASE1(:,npulse:n) = sin(walkphase1(:,npulse:n));
AVGX1 = mean(XPHASE1);
AVGY1 = mean(YPHASE1);
MAGXY1 = ( (AVGX1.^2) + (AVGY1.^2) ).^0.5;

plot(timestep,MAGXY1) % plot of this new average phase against time.
axis([0.9 2.0 0 1])
xlabel('Time (s)')
ylabel('Magnetization (A/m)')
%
% end
%
% plot to see how field gradient effects the magnitude of the spin echo.
Lecho = linspace(0,4,41);
Techo = exp(-1*(gamma)^2 *([Lecho]).^2* ((noff-Ninitial)*dt)^2 * 0.282 * (((non-noff)*dt) - ((noff-Ninitial)*dt)/3 ));
plot(Lecho,Techo)
TechoSIM = [1,0.81,0.38,0.16,0.03,0.005,0];
Lecho1 = [0,0.5,1,1.5,2,2.5,3];
hold on 
plot(Lecho1,TechoSIM,'r*')
xlabel('Field gradient (g/cm)')
ylabel('Echo Amplitude')

    %% constants.
Hbar = 6.63e-34/(2*pi);
Gamma = 42.58; %MHz/T
Bzero = 9.4; %T
K = 1.38e-23; %J/K
Temperature = 310; %K
DeltaE = Hbar*Gamma*Bzero; %J
SpinR = exp(-DeltaE/(K*Temperature));
% SpinR = (DeltaE/(K*Temperature))*0.5*N;

