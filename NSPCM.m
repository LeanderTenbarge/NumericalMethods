
% Leander Tenbarge
% Script to determine the Neumann solution to the one-dimensional, one-phase, Stefan problem:
clear all

% Actual Calculations: %%%
% Initial parameters:
To = 300;        % Boundary temperature @ x=0 [K]
Tm = 273;        % Melting temperature [K]
Cl = 2100;       % Specific heat capacity [J/kg*K]
L = 334000;      % Latent Heat of Fusion [J/kg]
alphal = 1.4e-7; % Thermal diffusivity [m^2/s]
LengthTime = 1e+7;
nt = 100;

% Calculating the Stefan Number:
Stl = (Cl*(To - Tm))/L; 

% Determining lambda (similarity parameter)
f = @(lambda) lambda*exp(lambda^2)*erf(lambda)*sqrt(pi) - Stl;
l = fzero(f,0.1);

% Interface position function
s = @(t) 2*l*sqrt(alphal*t);

% Temperature field function for x <= interface
Tfield = @(x,t) To - (To - Tm) * (erf(x./(2*sqrt(alphal*t)))/erf(l));

% Spatial and temporal domain for the first plot
t = linspace(1,LengthTime,nt);   % time steps
x_max = 1;                   % max x for plotting
xdomain = linspace(0, x_max, 500);  % spatial points for full domain

% Generating the interface position for the second plot:
interfacePlotPosition = zeros(size(t));
for i = 1:length(t)
    interfacePlotPosition(i) = s(t(i));
end

interfacePlotVelocity = diff(interfacePlotPosition);

% Plotting the Animation %%%%%

% Trick for Plotting Animations:
figure
hold on

% Second subplot (interface vs time, static)
subplot(3,1,2)       % make this subplot active
plot(interfacePlotPosition, t, 'r', 'LineWidth', 1.5)  % static plot
grid on
xlabel('x [m]')
ylabel('t [sec]')
title('Interface Position vs. Time')


% Third subplot (Interface Velocity (Ds/dt) vs. time)
subplot(3,1,3);
plot( t(1:99),interfacePlotVelocity, 'g', 'LineWidth', 1.5)  
ylabel('V [m/s]')
xlabel('t [sec]')
title('Interface Velocity vs. Time')

% First subplot animation (temperature profile)
fig1 = subplot(3,1,1);
xlabel(fig1,'x [m]')
ylabel(fig1,'Temperature [K]')
title(fig1,'Temperature Profile [K/M]')
ylim(fig1,[Tm-0.5 To+.5])
grid on
hold on

% Create plots in this axes
interface_plot = plot(fig1, 0, Tm, 'ro', 'MarkerFaceColor','r'); % moving interface
Tplot = plot(fig1, NaN, NaN, 'b', 'LineWidth',2); % initial temperature distribution

for i = 1:length(t)
    time = t(i);
    
    % Current interface
    interfacePosition = s(time);
    
    % Vectorized temperature over the full domain
    temperatureDistribution = Tm*ones(size(xdomain));
    insideInterface = xdomain <= interfacePosition;
    temperatureDistribution(insideInterface) = Tfield(xdomain(insideInterface), time);
    
    % Update plots
    set(Tplot, 'XData', xdomain, 'YData', temperatureDistribution)
    set(interface_plot, 'XData', interfacePosition, 'YData', Tm)
    
    drawnow
    pause(0.05)
end