% to copute the stable and unstable dairations of \Gamma^+

clear

%% the stable diractions by the Poincaré map ( about 1 mint)
% par_aut = [0.2 0.2 5.7];
% T_per   = 5.881088436121394;
%
% gamma    =   [-5.242526776906579
%               -5.242527133727762
%                0.017632169555805];
% gammafix =  [gamma(1), gamma(3)];
% wdot        = @(t,var,T)(T*Rossfun_aut(var,par_aut));
% M           = @(init)IVP_LE(@(t,var)wdot(t,var,T_per),...
%               [0 1],init,1e-5);
% J           = MyJacobian(M,gamma);
% eigvalue    = eig(J);
% [eigvect,~] = eig(J);
% un1         = eigvect(:,1);
% ce1         = eigvect(:,2);
% st1         = eigvect(:,3);
% st          = [st1(1),st1(3)];
% un          = [un1(1),un1(3)];
% ce          = [ce1(1),ce1(3)];
%
% stref1       = gammafix + (0.0025*(st));
% stref2       = gammafix - (0.0025*(st));
% unref1       = gammafix + (0.6*(un));
% unref2       = gammafix - (0.6*(un));
%% The periodic orbit \Gamma

% initcond_per =  gamma;
% ivpinit      = ode45(@(t,var)wdot(t,var,T_per),[0 1],initcond_per);
% temp         = @(s)(deval(ivpinit,s));
% si           = linspace(0,1,200);
% solinit      = bvpinit(si,temp,T_per);
% BVPper       = per_Ros(par_aut,solinit);


%%

%load('detaold\st_un_ce.mat')


%% Shooting on varing r 
% this section takes between 10 sec up to 3 mints depends on rstep


% system parameters
b             =  0.2;
c             =  5.7;
amin          = -0.2;
amax          =  0.2;

% fixed point of the return map
gammafix = [-5.242526776906579   0.017632169555805];

% stable eigenvector of gamma.
st = [0.097459720341894   0.995182696264705];

% initial conditon of the shooting, close to the past equlibrium.
initcond_conn = [-0.007008925432331;...
    -0.035044627161653;...
    0.035044627161653]; %

% To preduce Fig~5 please use the respective Limit vector to
% effectivally find the last intersection time of the pullback attractor
% with \Sigma and the respective step size of r

% 1) for T = 125, Limits = [0 5], and rstep = 0.001.
% 2) for T = 135, Limits = [0 5], and rstep = 0.001.
% 3) for T = 145, Limits = [0,5], and rstep = 0.0005.
% 4) for T = 155, Limits = [3 7], and rstep = 0.0001.

T             = 125;
rstep         = 0.001;
tspan         = [-30,T];
h             = 1.5e-2;
% nn            = 1;
Limits        = [0 5];
rr            = 0.9:rstep:1;
eta           = NaN(1,length(rr));
T0            = NaN(1,length(rr));
for nn = 1: length(rr)
    r             =  rr(nn);
    par           = [b,c,r,amin,amax];
    odefun_nonaut = @(t,var)(Rossfun_nonaut(t,var,par));
    
    [~,t,var]     = IVP_Ross(odefun_nonaut,tspan,initcond_conn ,h);
    
    % interpolation
    x       = @(tt)interp1(t,var(1,:),tt,'spline');
    y       = @(tt)interp1(t,var(2,:),tt,'spline');
    z       = @(tt)interp1(t,var(3,:),tt,'spline');
    
    % the last intersection with the Poincaré section
    T0(nn) = LastIntersectonTime(t,var,Limits);
    
    eta(nn) = (([x(T0(nn)),z(T0(nn))] - gammafix)*st')/(st*st');
    X(nn,:) = var(1,:);
    Y(nn,:) = var(2,:);
    Z(nn,:) = var(3,:);
    disp(nn)
end

% plotting eta

FIG_eta = figure(1);
plot(rr,eta,'b-',rr,zeros(1,length(rr)),'r-','LineWidth',3);
xlim([0.9,1]);
ylim([-0.2 0.4])
set(0,'defaulttextInterpreter','latex')
xlabel('$r$')
ylabel('$\eta$','Rotation',0)
set(gca,'FontSize', 22);
set(gcf,'Position',[0.3e3 0.3e3 1e3 0.35e3])
box on

% saving pdf
set(FIG_eta,'Units','Inches');
pos1 = get(FIG_eta,'Position');
set(FIG_eta,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
    [pos1(3), pos1(4)])
print(FIG_eta,'Ross_r_vs_del_T155','-dpdf','-r0')

%% Plotting the connection at certain rates Fig~4 in the paper

load('data/eta_movie_155.mat')  % calling stored data

% two values of r for two diffrent connections 
r             = 0.9202212159423;
% r             = 0.995651959127;


tspan         = [-30,200];
h             = 1.5e-2;

Limits        = [0 10];

par           = [0.2,5.7,r,-0.2,0.2];
odefun_nonaut = @(t,var)(Rossfun_nonaut(t,var,par));

[~,t,var]     = IVP_Ross(odefun_nonaut,tspan,initcond_conn ,h);

x       = @(tt)interp1(t,var(1,:),tt,'spline');
y       = @(tt)interp1(t,var(2,:),tt,'spline');
z       = @(tt)interp1(t,var(3,:),tt,'spline');

% T0 = LastIntersectonTime(t,var,Limits);
% eta = (([x(T0),z(T0)] - gammafix)*st')/(st*st');

% plotting

FIG1=figure(2);
set(0,'defaulttextInterpreter','latex') %latex axis labels
hold on
plot3(var(1,:),var(2,:),var(3,:),'-b','LineWidth',2)
plot3(initcond_conn(1),initcond_conn(2),...
    initcond_conn(3),'r.','MarkerSize',20)
plot3(BVPper.y(1,:),BVPper.y(2,:),BVPper.y(3,:),'r-.','LineWidth',3)
xlabel('$x(t)$')
ylabel('$y(t)$')
zlabel('$z(t)$','Rotation',0)
box on
set(gca,'FontSize',20)
grid on
view([-40 20])
xlim([-15 15]);
ylim([-15 15]);
zlim([-1 25]);

% saving pdf
set(gcf,'Position',[0.04e3 0.04e3 0.68e3 0.6e3])
set(FIG1,'Units','Inches');
pos1 = get(FIG1,'Position');
set(FIG1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
    [pos1(3), pos1(4)])
print(FIG1,'Ross_conn_3d_r09202212159423','-dpdf','-r0')

FIG2=figure(3);
plot(t,var(2,:),'-b','LineWidth',2)
hold on
xlabel('$t$')
ylabel('$x(t)$','Rotation',0)
box on
set(gca,'FontSize',20)

xlim([-30 200])
ylim([-15 15])
grid on

% saving pdf
set(gcf,'Position',[0.5e3 0.3e3 0.68e3 0.3e3])
set(FIG2,'Units','Inches');
pos2 = get(FIG2,'Position');
set(FIG2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
    [pos2(3), pos2(4)])
print(FIG2,'Ross_conn_2d_r0995651959127','-dpdf','-r0')

%% Movie zeros of eta gives connection
load('data/eta_movie_155.mat')  % calling stored data 
figure
set(gcf,'Color','w')
V = VideoWriter('conn_vs_r155_20','MPEG-4');
V.FrameRate = 20;
open(V);

for nn = 1: length(rr)
    % interpolation
    x_conn  = @(tt)interp1(t,X(nn,:),tt,'spline');
    y_conn  = @(tt)interp1(t,Y(nn,:),tt,'spline');
    z_conn  = @(tt)interp1(t,Z(nn,:),tt,'spline');
    sss     = linspace(-10,T0(nn),5000);
    
    %subplot1
    subplot(6,5,[2:4,7:9,13:14])
    
    % 1) the periodic orbit;
    plot3(BVPper.y(1,:),BVPper.y(2,:),BVPper.y(3,:),'-r','LineWidth',2);
    hold on
    
    % 2) the past fixed point;
    plot3(initcond_conn(1),initcond_conn(2),initcond_conn(3),'r.','MarkerSize',20);
    
    % 3) the section;
    plot3(-8:1:-2,-8:1:-2,zeros(1,7),'-k','LineWidth',3);
    
    % 4) the connection;
    plot3(x_conn(sss),y_conn(sss),z_conn(sss),'b-','LineWidth',1);
    
    % 5) the end point of the conniction;
    plot3(x_conn(T0(nn)),y_conn(T0(nn)),z_conn(T0(nn)),'b.','MarkerSize',20);
    
    % 6) the fixed point of the return map (gammafix)
    plot3(gamma(1),gamma(2),gamma(3),'r.','MarkerSize',20);
    
    % 7) Figure settings.
    set(0,'defaulttextInterpreter','latex')
    axis([-10 12 -12 10 -0.5 30])
    yticks([-10  0  10]);
    xticks([-10 -5 0 5 10]);
    xlabel('$x$')
    ylabel('$y$','Rotation',0)
    set(gca,'FontSize', 22);
    box on
    grid on
    hold off
    view([10 30])
    % subplot2
    subplot(6,5,21:30)
    
    % 1) eta and its zeros
    plot(rr,eta,'b-',rr,zeros(1,length(rr)),'r-','LineWidth',3);
    hold on
    
    % 2) moving horizontal line
    plot(rr(nn)*ones(1,10),linspace(-0.2,0.4,10),'--','Color',...
        [0.65 0.65 0.65],'LineWidth',2)
    
    % 3) Figure settings
    xlim([0.9,1]);
    ylim([-0.2 0.4])
    set(0,'defaulttextInterpreter','latex')
    xlabel('$r$')
    ylabel('$\eta$','Rotation',0)
    set(gca,'FontSize', 22);
    set(gcf,'Position',[0.01e3 0.045e3 0.68e3 0.74e3])
    box on
    hold off
    
    drawnow ()
    
    
    F(nn) = getframe(gcf);
    writeVideo(V,F(nn));
end
close(V);
winopen('conn_vs_r155_20.mp4')

%% Functions

% ode of nonautonomous Rossler system
function [ dvar] = Rossfun_nonaut(t,var,par)

% ODE of Rossler system (Rossler 1976)
%   Inputs:
%           var: the variables x,y,z,a
%           par: the parameters b,c, amin and a max
%   Output:
%           dvar: the drivetive of var after one time step


x   = var(1);
y   = var(2);
z   = var(3);

b    = par(1);
c    = par(2);
r    = par(3);
amin = par(4);
amax = par(5);
Max  = amax - amin;

a      = @(t)(((Max/2).*((tanh((Max.*r*t/2))+1)))+amin);

dx = - y - z;
dy =   x + a(t).*y;
dz =   b + z*(x-c);



dvar=[dx;dy;dz];

end

% a function to find the last intersection time between the pullback
% attrctor and Poincare section
function [T0] = LastIntersectonTime(time,var,Limit)

x       = @(tt)interp1(time,var(1,:),tt,'spline');
y       = @(tt)interp1(time,var(2,:),tt,'spline');
sect    = @(tt)(y(tt)- x(tt)); % section x=y or x-y = 0
cond    = 0;
T0      = max(time);
while or(or(cond >= 0 , T0 >= (max(time) - Limit(1))),...
        (T0<=max(time)-Limit(2)))
    Tinit   = max(time) - (6*rand);
    T0      = NR(sect,Tinit);
    cond    = x(T0);
end
end


% newton raphson iteration
function [p,ind,conv] = NR(f,p0)
%NR is a function to find the zeros of a function by using Newton
%Raphson ittration.
%   Input:
%           1) a function f: which is the function that you want to solve
%           2) an intial guss p0
%   Output
%           1) p: the zero of the function f (i.e: f(p)=0)

PP=repmat(p0,1,20);
tol=1e-6; ind=0;
conv=0;
while and(conv<1,ind<20)
    ind=ind+1;
    df=MyJacobian(f,PP(:,ind));
    PP(:,ind+1)=PP(:,ind)-df\f(PP(:,ind));
    diff=norm(PP(:,ind+1)-PP(:,ind));
    root=norm(f(PP(:,ind)));
    if and(diff>=tol,and(root>=tol,ind<200))
        conv=0;
    else
        conv=1;
    end
    p=PP(:,ind);
end
end

% Jacobian of a n dimentional function 
function df = MyJacobian(f,x)
% to coumpute a jacobian matrix for a function f:R^n --> R^n
h=1e-12;
n=length(f(x)); m=length(x);
df = NaN(n,m);
if isfinite(sum(x))
    % F=repmat(f(x),1,m);
    for j=1:m
        F1=f(x);
        x(j)=x(j)+h;
        F2=f(x);
        x(j)=x(j)-h;
        df(:,j)=(F2-F1)./h;
    end
else
    df = NaN(n,m);
end
end

% IVP solver 

function [xend,t,xt]=IVP_Ross(f,tspan,x0,h)

%[t,xt] = IVP(f,x0,tspan,h)
%   Inputs:
%       f: function defining the right-hand side of the ODE. f should accept two arguments: t (a number) and x (an n-dimensional vector). The function f should return an n-dimensional vector y (the time derivative). Typical calling sequence: y=f(t,x), returning the value of f at time t in position x.
%       x0: initial value where integration starts from (n-dimensional vector)
%       tspan: Starting time and end time for integration. Integration has to run from time t =tspan(1) to time t =tspan(2).
%       N: number of steps for integration. The integration stepsize h=(tspan(2)-tspan(1))/N should be small.
%   Outputs:
%       t: vector of times at which intermediate values have been computed (this should have N + 1 entries).
%       xt: intermediate values (n  (N + 1)-array). xt(:,k) should be the solution at t(k).
N=fix(abs(tspan(2)-tspan(1))/abs(h));
if tspan(1)<tspan(2)
    h=abs(h);
else
    h=-1*abs(h);
end
t = tspan(1):h:tspan(2);
xt(:,1) = x0;

for i=1:N
    y1=f(t(i),xt(:,i));
    y2=f(t(i)+h/2,xt(:,i)+y1*h/2);
    y3=f(t(i)+h/2,xt(:,i)+y2*h/2);
    y4=f(t(i)+h,xt(:,i)+y3*h);
    D=h*(y1+2*y2+2*y3+y4)/6;
    xt(:,i+1)=xt(:,i)+D;
    xend=xt(:,i+1);
end
% xt = xt(:,1:end-1);
end
