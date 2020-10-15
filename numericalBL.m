%% ------------------------------------------------------------------------
% Numerical Solution of Similar Boundary Layer Equations.
% 
% 
% author./ dekeract01
% 
% -------------------------------------------------------------------------
clear all; close all; clc; warning off
tic
colorhandles; global cyan brown blue green red 
%% The Blasius Boundary Layer on a Flat Plate at Zero Pressure Gradient.
% Blasius Equation
% f^{111} + f f^{11} = 0.
% boundary conditions@: f^{1}(\eta = 0) = 0
%                       f(\eta = 0)     = 0
%                       f^{1}(\eta \rightarrow \infty) = 1
% -------------------------------------------------------------------------
tic
initial_conditions = [0,0,0.33];
vwrange=-1:0.1:1;
birange=[-1 3];
conval=1;
tol=1e-10;
beta=0;
position=2;
changeconditions=[vwrange',zeros(size(vwrange,2),1),...
    zeros(size(vwrange,2),1)];

% Solve Blasius Equation Numerically
[xs,ys] = blasius(initial_conditions,beta,'Blasius');
% Plot results of the Blasius Equation
figure(1); hold on;
plot(ys(:,1),xs,'-','Linewidth',2,'Color',brown)
plot(ys(:,2),xs,'-','Linewidth',2,'Color',red)
plot(ys(:,3),xs,'-','Linewidth',2,'Color',blue)
title('Numerical Solution of the Blasius Equation')
figformat


% Effect of Wall Transpiration - bisection method
for i=1:size(vwrange,2)
    [conconditions(i,:)] = bisection(initial_conditions,vwrange(i),...
        birange,conval,tol,beta,'Blasius',position);
end

% Effect of Wall Transpiration -fsolve test
for i=1:size(vwrange,2)
    [concfsolve(i,:),~,~] = bvpsolve(initial_conditions,beta,'Blasius',...
        changeconditions(i,:));
end

figure(2); hold on;
plot(concfsolve(:,3),concfsolve(:,1),'-','Linewidth',2,'Color',blue)
title('Parametric Study of Wall Transpiration of Blasius Equation vs v_w')
xlabel('f(0) = v_w ');
ylabel('f^{''}(0) = \tau_w');
grid on
grid minor
xlim([0 1.5]);
%xticks(0:0.5:2);
ylim([-1 1.5]);
% Plot results

%% Falkner and Skan
% Blasius Equation
% f^{111} + f f^{11} + \beta(1- (f^{11})^{2} ) = 0.
% boundary conditions@: f^{1}(\eta = 0) = 0
%                       f(\eta = 0)     = 0
%                       f^{1}(\eta \rightarrow \infty) = 1
% -------------------------------------------------------------------------
beta=0.05;
% Solve the Falkner-Skan Equation Numerically
[xb,yb] = blasius(initial_conditions,beta,'FlankerandSkan');
figure(3); hold on;
plot(yb(:,1),xb,'-','Linewidth',2,'Color',brown)
plot(yb(:,2),xb,'-','Linewidth',2,'Color',red)
plot(yb(:,3),xb,'-','Linewidth',2,'Color',blue)
title('Numerical Solution of the Falkner-Skan Equation, \beta=0.1')
figformat

% Parametric study: Pressure Gradient
betap=-0.2:0.1:0.5;
birange=[-1 3];
conval=1;
tol=1e-10;
vwbrange=-1:0.1:1;

for i=1:size(betap,2)
    if betap(i)<=-0.2
        tol=1e-5;
    elseif betap(i)>=0.7
        tol=1e-2;
    end
    [betaconditions(i,:)] = bisection(initial_conditions,0,birange,...
        conval,tol,betap(i),'FlankerandSkan',position);
end

betanew=0.05;
% Parametric study: Wall Transpiration
for i=1:size(vwbrange,2)
    [betatau_conconditions(i,:)] = bisection(initial_conditions,...
        vwbrange(i),birange,conval,tol,betanew,'FlankerandSkan',position);
end

% Displacement thickness \delta^{*}
for i=1:size(betatau_conconditions,1)
    deltastar(i,1)=betatau_conconditions(i,1);
    [deltastar(i,2)] = displcaementthickness(betatau_conconditions(i,:),beta,'FlankerandSkan');
end

figure(4)
plot(deltastar(:,1), deltastar(:,2), 'LineWidth', 2, 'color', blue);
title('Change in Displacement thickness \delta^{*} at differnt intial conditions f(0)=v_w', 'FontSize', 14);
xlabel('v_w', 'FontSize', 20);
ylabel('\delta^{*}', 'FontSize', 20);
axis tight
grid on
grid minor 

%% Total Skin friction coefficient
rho=1.26;
u=0.5;

Cf=betaconditions(:,3)*0.5*(1/(rho*u^2));


figure(5)
plot(Cf, -betap, 'LineWidth', 2, 'color', 'black');
title(' Skin friction coeff at against \beta', 'FontSize', 14);
xlabel('C_{f}', 'FontSize', 20);
ylabel('-\beta', 'FontSize', 20);
axis tight
grid on
grid minor 

U_inf = 0.5;
L = 10;
mu = 1.789E-5;
rho = 1.225;
nu = mu/rho;
A = sqrt(nu/U_inf);
h = 0.01;

eta=xb;
% Velocity Profile and the BL thickness distribution
figure(6)
hold all
for i = 1:length(xb)
delta(i) = 5*sqrt(xb(i))*A;
end
plot(xb, delta, 'LineWidth', 2, 'color', 'black');
position = [1 5 8];
for j = 1:length(position)
y{j} = eta*sqrt(position(j))*A;
plot(yb(:,2)+position(j), y{j}, 'LineWidth', 2)
end
Legend2 = {'\delta (x)', 'Velocity Profile at x = 1', 'Velocity Profile at x = 5', 'Velocity Profile at x = 8'};
title('Velocity profiles at x = 1, 5 and 8 plus BL thickness distribution along the plate', 'FontSize', 14);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
legend(Legend2, 'FontSize', 14);
ylim([0, 2*max(y{j})])
grid on




%% Illington and Stewartson
beta=0.1;
initial_conditions = [0,0,0.5,0,0.5];
%[xi,yi] = blasius(initial_conditions,beta,'Stewartson');
changeval=[0.5,0,0,0.5,0];
[~,xi,yi] = bvpsolve(initial_conditions,beta,'Stewartson',changeval);
posnew=5;
figure(7); hold on;
plot(yi(:,1),xi,'-','Linewidth',2,'Color',brown)
plot(yi(:,2),xi,'-','Linewidth',2,'Color',red)
plot(yi(:,3),xi,'-','Linewidth',2,'Color',blue)
plot(yi(:,5),xi,'-','Linewidth',2,'Color',green)
plot(yi(:,4),xi,'-','Linewidth',2,'Color',cyan)
title('Numerical Solution of the Compressible Boundary Layer Equation, \beta=0.1, v_w=0.5, S_w = 0.5')
figformatlast


betai=-0.1:0.05:0.1;
sw=-1:0.1:1;
constsw=0.8;
conval=0;
v_w=0.5;

changebeta=[v_w,0,0,constsw,0];
changesw=[v_w*ones(size(sw,2),1),zeros(size(sw,2),1),zeros(size(sw,2),1),sw',zeros(size(sw,2))];

% Parametric Studt with varying \beta
for i=1:size(betai,2)
    [stebetafsolve(i,:),~,~] = bvpsolve(initial_conditions,betai(:,i),'Stewartson',changebeta);
end

% Parametric Study with varying S_w 
for i=1:size(sw,2)
    [steswfsolve(i,:),~,~] = bvpsolve(initial_conditions,beta,'Stewartson',changesw(i,:));
end

figure(8)
plot(steswfsolve(:,4), steswfsolve(:,3), 'LineWidth', 2, 'color', cyan);
title(' Variation of \tau_w with a change in the S_w values', 'FontSize', 14);
xlabel('s_w', 'FontSize', 20);
ylabel('\tau_w', 'FontSize', 20);
axis tight
grid on
grid minor 


toc

%% Functions
% -------------------------------------------------------------------------
% Function. [x,y] = blasius(initial_conditions,beta,whichone)
%           The following function takes a variable and finds it's initial
%           values regardless of the guess and returns the x and y
%           matrices. 
%           # initial_conditions: guess of initial conditions. 
%           # beta: only applicable for Flakner and Skan 
%           # whichone: choose the boundary layer case. 
%
% author./ t.ala. 2019
% university of liverpool.
% -------------------------------------------------------------------------
function [x,y] = blasius(initial_conditions,beta,whichone)
% Use fsolve to ensure the boundary function is zero. The result is the
% unknown initial condition.
opt = optimset('Display','off','TolFun',1E-20);
F = fsolve(@(F) eval_boundary(F,beta,whichone),initial_conditions,opt);
if strcmpi(whichone,'Stewartson')
    %initial_conditions(5)=0;
    %F = fsolve(@(F) eval_boundary(F,beta,whichone),initial_conditions,opt);
    %F(1)=0;
    %F(2)=0;
end
% Solve the ODE-IVP with the converged initial condition
[x,y] = solve_ode(F,beta,whichone);
end

function [F,x,y] = bvpsolve(initial_conditions,beta,whichone,...
    changeconditions)
% Use fsolve to ensure the boundary function is zero. The result is the
% unknown initial condition.
opt = optimset('Display','off','TolFun',1E-20);
F = fsolve(@(F) solveini(F,beta,whichone,changeconditions)...
    ,initial_conditions,opt);
if strcmpi(whichone,'Stewartson')
    %initial_conditions(5)=0;
    %F = fsolve(@(F) eval_boundary(F,beta,whichone),initial_conditions,opt);
    %F(1)=0;
    %F(2)=0;
end
% Solve the ODE-IVP with the converged initial condition
[x,y] = solve_ode(F,beta,whichone);
end

function [delta] = displcaementthickness(initial_conditions,beta,whichone)
% Use fsolve to ensure the boundary function is zero. The result is the
% unknown initial condition.
opt = optimset('Display','off','TolFun',1E-20);
F = initial_conditions;
% Solve the ODE-IVP with the converged initial condition
[x,y] = solve_ode(F,beta,whichone);
delta = x(end)-y(end,1) - (x(1)-y(1,1));
end

function [conconditions] = bisection(initial_conditions,...
    vwrange,birange,conval,toler,beta,whichone,position)
% Use bisection technique to check for boundary conditions at various
% points
iter=0;
for i=1:size(vwrange,2)
    initial_conditions(1)=vwrange(i);
    F = initial_conditions;
    [~,y] = solve_ode(F,beta,whichone);
    median=1;
    while abs(1-y(end,position))>toler
        iter=iter+1;
        median=0.5*(birange(2)+birange(1));
        if strcmpi(whichone,'Blasius') || strcmpi(whichone,'FlankerandSkan')
            F = [initial_conditions(1) initial_conditions(2) birange(1)];
            F1 = [initial_conditions(1) initial_conditions(2) birange(2)];
            Fm = [initial_conditions(1) initial_conditions(2) median];
        else
            F = [initial_conditions(1) initial_conditions(2) ...
                initial_conditions(3) initial_conditions(4) birange(1)];
            F1 = [initial_conditions(1) initial_conditions(2) ...
                initial_conditions(3) initial_conditions(4) birange(2) ];
            Fm = [initial_conditions(1) initial_conditions(2) ...
                initial_conditions(3) initial_conditions(4) median];
        end
        
        [~,y] = solve_ode(F,beta,whichone);
        [~,y1] = solve_ode(F1,beta,whichone);
        [~,ym] = solve_ode(Fm,beta,whichone);
        if abs((conval-y(end,position))) < toler
            y=y;
            m=birange(1);
            conconditions(i,:)=[initial_conditions(1) initial_conditions(2)...
                initial_conditions(3) initial_conditions(4) m];
            return
        elseif abs((conval-y1(end,position))) < toler
            y=y1;
            m=birange(2);
            conconditions(i,:)=[initial_conditions(1) initial_conditions(2)...
                initial_conditions(3) initial_conditions(4) m];
            return
        elseif abs((conval-ym(end,position))) <toler
            y=ym;
            m=median;
            conconditions(i,:)=Fm;
            return
        end
        rootsign1=conval-y(end,position);
        rootsign2=conval-y1(end,position);
        rootsign3=conval-ym(end,position);
        rootsign1=rootsign1/abs(rootsign1);
        rootsign2=rootsign2/abs(rootsign2);
        rootsign3=rootsign3/abs(rootsign3);
        
        if rootsign3 == rootsign1
            birange(1)= median;
            y=y;
        elseif rootsign3 == rootsign2
            birange(2)= median;
            y=y;
        end
        
    end
end
end
function [x,y] = solve_ode(F,beta,whichone)
% Solve the ODE-IVP with initial condition F on [0 100] (arbitrary upper
% bound)
if strcmpi(whichone,'Blasius')
    [x,y] = ode45(@(x,y) [y(2); y(3); -y(1)*y(3)],[0 10],F); %solve BVP
elseif strcmpi(whichone,'FlankerandSkan')
    [x,y] = ode45(@(x,y) [y(2); y(3); -y(1)*y(3)-(beta*(1-y(2)^2))],...
        [0 10],F); %solve BVP
elseif strcmpi(whichone,'Stewartson')
    [x,y] = ode45(@(x,y) [y(2); y(3); -y(1)*y(3)+(beta*(((y(2)*y(2))...
        -1-y(4)))); y(5); -y(1)*y(5)],[0 10],F); %solve BVP
end
end

function [g] = eval_boundary(F,beta,whichone)
% Get the solution to the ODE with inital condition F
[x,y] = solve_ode(F,beta,whichone);
% Get the function values (for BCs) at the starting/end points
f_start = y(1,1); %f(0) = 0 || defined
df_start = y(1,2); %f'(0) = 0
df_end = y(end,2); %f'(inf) - 1 = 0


% Evaluate the boundary function
g = [f_start - 0.5
    df_start
    df_end - 1];

if strcmpi(whichone,'Stewartson')
    f_start = y(1,1); %f(0) = 0 || defined
    df_start = y(1,2); %f'(0) = 0
    df_end = y(end,2); %f'(inf) - 1 = 0
    s_start = y(1,4); % s(0) = 0 || defined
    ds_end = y(end,4); % s(inf) = 0
    % Evaluate the boundary function
    g = [f_start
        df_start
        df_end - 1
        s_start
        ds_end];
end

end
function figformat
% This function simply formats the figure
xlabel('\it{f, f^{(1)}, f^{(2)}}');
ylabel('\eta');
xlim([-0.5 2]);
%xticks(0:0.5:2);
ylim([0 10]);
h = legend('\it{f}','\it{f}^{ (1)}','\it{f^{ (2)}}',...
    'Location','NorthEast');
set(h,'FontSize',14);
axis square
set(gca, 'box', 'on'); %creates a box border
set(gca,'FontWeight','normal','linewidth',1.05);...
    fontname = 'Arial'; %the next few lines format the plot
set(0,'defaultaxesfontname',fontname);...
    set(0,'defaulttextfontname',fontname);
fontsize = 16;
set(0,'defaultaxesfontsize',fontsize);...
    set(0,'defaulttextfontsize',fontsize);
grid on
grid minor
end
function figformatlast
% This function simply formats the figure
xlabel('\it{f, f^{ (1)}, f^{ (2)},S,S^{ (1)}}');
ylabel('\eta');
xlim([-1 2]);
%xticks(0:0.5:2);
ylim([0 10]);
h = legend('\it{f}','\it{f}^{ (1)}','\it{f^{ (2)}}',...
    '\it{S}','\it{S}^{ (1)}','Location','NorthEast');
set(h,'FontSize',18);
axis square
set(gca, 'box', 'on'); %creates a box border
set(gca,'FontWeight','normal','linewidth',1.05);...
    fontname = 'Arial'; %the next few lines format the plot
set(0,'defaultaxesfontname',fontname);...
    set(0,'defaulttextfontname',fontname);
fontsize = 20;
set(0,'defaultaxesfontsize',fontsize);...
    set(0,'defaulttextfontsize',fontsize);
grid on
grid minor
end
% -------------------------------------------------------------------------
%
%
%
%
% -------------------------------------------------------------------------
function g = solveini(F,beta,whichone,changeconditions)
% Get the solution to the ODE with inital condition F
[x,y] = solve_ode(F,beta,whichone);
% Get the function values (for BCs) at the starting/end points
f_start = y(1,1); %f(0) = 0 || defined
df_start = y(1,2); %f'(0) = 0
df_end = y(end,2); %f'(inf) - 1 = 0



% Evaluate the boundary function
g = [f_start + (-1*changeconditions(1))
    df_start
    df_end - 1];

if strcmpi(whichone,'Stewartson')
    f_start = y(1,1); %f(0) = 0 || defined
    df_start = y(1,2); %f'(0) = 0
    df_end = y(end,2); %f'(inf) - 1 = 0
    s_start = y(1,4); % s(0) = 0 || defined
    ds_end = y(end,4); % s(inf) = 0
    % Evaluate the boundary function
    g = [f_start + (-1*changeconditions(1))
        df_start
        df_end - 1
        s_start + (-1*changeconditions(4))
        ds_end];
end

end

function colorhandles
% figure colours
global cyan brown orange blue green red 
cyan        = [0.2 0.8 0.8];
brown       = [0.2 0 0];
orange      = [1 0.5 0];
blue        = [0 0.5 1];
green       = [0 0.6 0.3];
red         = [1 0.2 0.2];
end

% -------------------------------------------------------------------------
% the end ...
% -------------------------------------------------------------------------
