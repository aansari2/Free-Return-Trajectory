%%% Free-Return Trajectory Simulation
function FreeReturn
colormap(winter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global data G name
data = [
    % Mass   MinRadius(km)  InitVel(km/s)  Radius(km)   Object
    5.972*10^24   4414.6      0.0133        6378.1      % Earth
    7.342*10^22   406000      0.957184      1738.1      % Moon
    48*10^00      6655        10.88         500         % Satellite
    ];
n = size(data,1); % Number of Bodies
name =  {'Earth','Moon','Craft'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = 6.67408 * 10^-11; % m3 kg^-1 s^-2 (Gravitational Constant)
% F_vector = GMm/r^3*r_vector (Gravity Equation)
%%%%% Initial Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = [pi 0 -2.30885]'; % angular configuration
% Inital conditions
IC = [
    data(:,2)*1e3.*cos(config);   % X-Positions
    data(:,2)*1e3.*sin(config);   % Y-Positions
    -data(:,3)*1000.*sin(config); % X-Velocities
    data(:,3)*1000.*cos(config)]; % Y-Velocities
%%%%%%%%%%%%%%%(Ensure sum_Momentums = 0)%%%%%%%%%%%%%%%%
IC([0 1]*n+2) = IC([0 1]*n+2) + IC([0 1]*n+1);
IC([0 1]*n+3) = IC([0 1]*n+3) + IC([0 1]*n+1);
IC(3*n+1) = -sum(IC(3*n+2:end).*data(2:end,1)/data(1,1)); %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 60000;  % Number of intervals
dt = 10;    % time step length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0;                % Initial Time
pos = zeros(4*n,N+1); % Preallocate position
pos(:,1) = IC;        %
plot(0);clf
for k = 1:N
    %%%%% This is a psuedo energy-preserving Runge kutta scheme
    % Runge-Kutta 3 Stages
    k1 = dt*rhs(t, pos(:,k));                                   % Stage 1
    k2 = dt*rhs(t + 1/3*dt, pos(:,k) + 1/3*k1);                 % Stage 2
    k3 = dt*rhs(t + 5/6*dt, pos(:,k)+(-5/48*k1 + 15/16*k2));    % Stage 3
    
    pos(:,k+1) = pos(:,k) + (1/10*k1 + 1/2*k2 + 2/5*k3);        % Integrate position, velocity
    t = t + dt;                                                 % Integrate t
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Plot every 1000 iterations or at the end
    if mod(k,1000)==1 || k==N
        dims = 5e8*[-.1 1 -.1 1];
        X = linspace(dims(1),dims(2),128*2)';
        Y = linspace(dims(1),dims(2),128*2) ;
        U = X*0+Y*0;
        for i = 1:2
            R = sqrt((pos(i,k)-(Y*0+X)).^2 + (pos(i+n,k)-(Y+X*0)).^2);
            U = U - G*data(i,1)./R;
        end
        imagesc(X,Y,U'); caxis([-5e6 0]); hold on;
        for i = 1:n
            plot(pos(i,1:k+1),pos(i+n,1:k+1))
            hold on;
        end
        col = {'w';'k';'k'};
        for i = 1:n
            ck = exp(2i*pi*(0:0.05:1))*data(i,4)*1e3 + pos(i,k+1) + 1i*pos(i+n,k+1);
            fill(real(ck),imag(ck),[100,202,157]/255,'facealpha',.5)
            text(pos(i,k+1),pos(i+n,k+1),['  ' name{i}],'color',col{i})
        end
        title(num2str(t/3600/24,'t = %.2f days | Free Return Trajectory'))
        hold off
        axis equal;
        xlabel('x-axis (meters)');
        ylabel('y-axis (meters)');
        set(gca,'Ydir','normal');
        axis(dims)
        ylabel(colorbar,'Gravitaitonal Potential');
        drawnow;
    end
end
%%%% Find Acceleration as a function of position for n-bodies %%%%%%%%%%%%%%
    function dzdt = rhs(~,z)
        dzdt = zeros(4*n,1);
        dzdt(1:2*n) = z(2*n+1:end); % derivative of position is velocity
        
        %%%%%% Compute Accelerations of Each object %%%%%%%%%%%
        index = 1:n;
        mass = data(:,1);
        x = z(index);     % x - Position
        y = z(index+n);   % y - Position
        a_x = zeros(n,1); % x - acceleration preallocate
        a_y = zeros(n,1); % y - acceleration preallocate
        for j = 1:n
            % serpation distances
            r = ((x - x(j)).^2+(y - y(j)).^2).^.5;
            
            % x - acceleration
            a_x(j) = sum(G*mass(index~=j)./r(index~=j).^3.*(x(index~=j) - x(j)));
            
            % y - acceleration
            a_y(j) = sum(G*mass(index~=j)./r(index~=j).^3.*(y(index~=j) - y(j)));
        end
        dzdt(2*n+1:end) = [a_x;a_y]; % derivative of velocity is acceleration
    end
end
% code by u/ansariddle
