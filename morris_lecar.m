clear 
close all
clc

%% code font settings
%%%% Set "Arial" as the Default font
set(0,'defaultAxesFontSize',16);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultTextFontSize',16);
set(0,'defaultTextFontName','Arial');

set(0,'defaultUipanelFontName','Arial');
set(0,'defaultUicontrolFontName','Arial');
%%
Nt     = 50000;  % Num. of sample
dt     = 0.01;   % time step for numerical integration; unit : msec
time   = linspace(0, Nt-1, Nt) * dt; % time vector; unit : msec
%%%%% parameter settings
%%% typical parameter setting for Type I mode


C    =  5;
gL   =  2;
gK   =  8;
gCa  =  4;
VL   = -60;
VK   = -80;
VCa  =  120;
V1   = -1.2;
V2   =  18;
V3   =  12;
V4   =  17.4;
Iext =  39.8;
phi  =  1/15;



%%% typical parameter setting for Type II mode
% C    =  5;
% gL   =  2;
% gK   =  8;
% gCa  =  4.4;
% VL   = -60;
% VK   = -80;
% VCa  =  120;
% V1   = -1.2;
% V2   =  18;
% V3   =  2;
% V4   =  30;
% Iext =  100.5;
% phi  =  1/25; %unit: 1/msec 

X0     = [0, 0]; % initial value of state variables
                 % X0(1): membrane potential, v
                 % X0(2): recovery variable,  w
%%%%% parameter settings
%% Solve differential equation
X      = zeros(Nt, length(X0));
X(1,:) = X0;

for i = 2:Nt
    X_now  = X(i-1,:);
    %%%%% Numerical integral scheme with 4th order Runge Kutta method
    X(i,:) = runge_kutta(X_now, dt, @MorrisLecar, ...
                                    C, gL, gK, gCa,...
                                       VL, VK, VCa,...
                                       V1, V2, V3, V4,...
                                       Iext, phi);
end
%%
fig = figure(1);
% figure_setting(60, 40, fig);

sfh1 = subplot(2,1,1,'parent', fig);
plot(time, X(:,1), 'LineWidth', 3);
hold on
plot(time, X(:,2), 'LineWidth', 3);
hold off

xlabel('time (ms)')
ylabel('V, N')
lgnd = legend({'membrane potential \it V', 'recovery variable \it N'}, 'location', 'northeastoutside');
%%%%%%%
sfh2 = subplot(2,1,2,'parent', fig);
plot(X(:,1), X(:,2), 'r', 'LineWidth', 3);
xlabel('membrane potential \it V')
ylabel('recovery variable \it N')
title('phase space')
axis square
sfh2.Position = sfh2.Position - [0.1, 0, 0, 0];

% fname = [filepath, filesep, 'figures', filesep, 'ex1', filesep, 'result'];
% figure_save(fig, fname)


function dXdt = MorrisLecar(X, varargin)
    V    = X(1);
    N    = X(2);
    
    if length(varargin)==1    
        par  = varargin{1};
    else
        par  = varargin;
    end
    
    C    = par{1};
    gL   = par{2};
    gK   = par{3};
    gCa  = par{4};
    VL   = par{5};
    VK   = par{6};
    VCa  = par{7};
    V1   = par{8};
    V2   = par{9};
    V3   = par{10};
    V4   = par{11};
    Iext = par{12};
    phi  = par{13}; 

    Minf = Sigm(V, V1, V2);
    Ninf = Sigm(V, V3, V4);
    
    dVdt = 1/C * (- gL  * (V - VL) ...
                  - gCa * Minf * (V - VCa) ...
                  - gK  * N  * (V - VK) + Iext);
    dNdt =  Lambda(V, V3, V4, phi) * (Ninf - N);

    dXdt = [dVdt, dNdt];
end

function val = Sigm(V, V1, V2)
    %%%% sigmoid function
    val =  1 / (1 + exp(-2 * (V - V1)/V2));
    % This function can be also expressed as: val = 0.5 * (1 + tanh((V - V1)/V2)); 
end

function lambda = Lambda(V, V1, V2, phi)
    lambda = phi * cosh((V-V1)/(2*V2));
end

function X_next = runge_kutta(X_now, dt, func, varargin)
    k1     = func(X_now, varargin);
    
    X_k2   = X_now + (dt/2) * k1;
    k2     = func(X_k2, varargin);
    
    X_k3   = X_now + (dt/2) * k2;
    k3     = func(X_k3, varargin);
    
    X_k4   = X_now + dt * k3;
    k4     = func(X_k4, varargin);

    X_next = X_now + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end