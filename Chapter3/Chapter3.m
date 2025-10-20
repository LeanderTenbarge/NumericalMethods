% Leander Tenbarge
% Chapter 3: Parabolic PDE's, Hoffman's Numerical Methods


clear all

% Input:
height = .040; % m
time = 5;    % sec
dy = 0.001;    % m
dt = 0.005;    % sec
nu = 0.000217; % m^2/s
Uo = 40;       % m/s 
r = nu*dt/dy^2;

% Generating spacial and temporal domain
y = 0:dy:height;                    % Spatial domain from 0 to height with step dx
t = 0:dt:time;                      % Temporal domain from 0 to 1 second with step dt
Uftcs = zeros(length(y),length(t)); % FTCS explicit solution
Udf = zeros(length(y),length(t));   % Dufort - Frankel method
Ucn = zeros(length(y),length(t));   % Crank-Nicholson implicit
Ul = zeros(length(y),length(t));    % Laasonen method



% FTCS explicit method:
for j = 2:length(t) % Looping over the temporal domain:
    Uftcs(1,:) = Uo;      % U(y = 0) = Uo
    Uftcs(end,:) = 0;     % U(y = h) = 0 

    for i = 2:length(y)-1 % Looping over the spacial domain
        Uftcs(i,j) = Uftcs(i,j-1) + ((nu*dt)/(dy^2))*(Uftcs(i+1,j-1) - 2*Uftcs(i,j-1) + Uftcs(i-1,j-1));
    end
end



% Dufort - Frankel method:
for j = 2:length(t) % Looping over the temporal domain:
    Udf(1,:) = Uo;      % U(y = 0) = Uo
    Udf(end,:) = 0;     % U(y = h) = 0 

    for i = 2:length(y)-1 % Looping over the spacial domain
        Udf(i,j+1) = ((1 - 2*r)*Udf(i,j-1) + 2*r*(Udf(i+1,j) + Udf(i-1,j)))/(1 + 2*r);
    end
end


% Crank - Nicholson:
% Contructing the A & B matricies:
A = zeros(length(y), length(y)); 
B = zeros(length(Ucn(:,1)), length(Ucn(1,:))-1); 
A(1,1) = 1; 
A(end,end) = 1;

for i = 2:length(Udf(:,1))-1
    for j = 2:length(Udf(:,1))-1
        if i ==j
            A(i,j-1) = -r;
            A(i,j) = 2*(1+r);
            A(i,j+1) = -r;
        end 
    end
end

B(1,:) = Uo;
B(end,:) = 0;

for j = 2:length(Ucn(1,:)) 
    for i = 2:length(Ucn(:,1)) - 1
        B(i,j) = r*Ucn(i-1,j-1) + 2 * (1 - r)*Ucn(i,j-1) + r*Ucn(i+1,j-1);
    end 
    Ucn(:,j) = A \ B(:,j);
end











% % FTCS explicit plotting:
% subplot(3,1,1)
% hold on 
% for i = [find(t == .18),find(t == 0.36),find(t == .54),find(t == .72),find(t == .90),find(t == 1.08),find(t == 5)]
%     plot(Uftcs(:,i),y,'^-')
% end
% title(sprintf('Velocity Profiles obtained by FTCS Explicit method, \\Delta t = %.3f s, \\Delta y = %.3f m', dt, dy))
% legend('t = .18', 't = 0.36', 't = 0.54', 't = 0.72', 't = .90', 't = 1.08' )
% xlabel('u (m/sec)')
% ylabel('y (m)')
% hold off
% 
% % Dufort - Frankel explicit plotting:
% subplot(3,1,2)
% hold on 
% for i = [find(t == .18),find(t == 0.36),find(t == .54),find(t == .72),find(t == .90),find(t == 1.08),find(t == 5)]
%     plot(Udf(:,i),y,'^-')
% end
% title(sprintf('Velocity Profiles obtained by Dufort - Frankel Explicit method, \\Delta t = %.3f s, \\Delta y = %.3f m', dt, dy))
% legend('t = .18', 't = 0.36', 't = 0.54', 't = 0.72', 't = .90', 't = 1.08' )
% xlabel('u (m/sec)')
% ylabel('y (m)')
% hold off
% 
% % Crank - Nicholson Plotting
% subplot(3,1,3)
% hold on 
% for i = [find(t == .18),find(t == 0.36),find(t == .54),find(t == .72),find(t == .90),find(t == 1.08),find(t == 5)]
%     plot(Ucn(:,i),y,'^-')
% end
% title(sprintf('Velocity Profiles obtained by Crank - Nicolson method, \\Delta t = %.3f s, \\Delta y = %.3f m', dt, dy))
% legend('t = .18', 't = 0.36', 't = 0.54', 't = 0.72', 't = .90', 't = 1.08' )
% xlabel('u (m/sec)')
% ylabel('y (m)')
% hold off


for i = [find(t == .18),find(t == 0.36),find(t == .54),find(t == .72),find(t == .90),find(t == 1.08),find(t == 5)]
    plot(Udf(:,i),y,'^-')
end
title(sprintf('Velocity Profiles obtained by Dufort - Frankel Explicit method, \\Delta t = %.3f s, \\Delta y = %.3f m', dt, dy))
legend('t = .18', 't = 0.36', 't = 0.54', 't = 0.72', 't = .90', 't = 1.08' )
xlabel('u (m/sec)')
ylabel('y (m)')
hold off
