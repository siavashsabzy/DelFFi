clear
clc
time0 = [2014 09 14 12 00 00];               % Satellite Start Mission
times = [2014 08 13 12 00 00];               % Sun last allingment
I     = [0.017,0,0;0,0.055,0;0,0,0.055];     % Satellite Moment Of Inertia
mu    = 398576.0576;                         % Gravitational Coe.
RE    = 6371;                                % in km
alt   = 200;                                 % altitude in km
a0    = (alt+RE);                            % semi-major axis, km
e0    = 0;                                   % eccentricity
i0    = 98.6;                                % inclination, deg
Raan0 = -15;                                 % right ascension of the ascending node, deg
aop0  = 0;                                   % argument of periapse, deg                           
ta0   = 0;                                   % true anomaly, deg
p0    = a0*(1-e0^2);                         % p0
nO    = sqrt(mu/(a0^3));                     % Orbit angular Rate
w0    = 5;                                   % degree per second
Bbias = [500,500,500]*10^-9;                 % magnetic field measurement bias
Bn    = [150,150,150]*10^-9;                 % magnetic field measurement noise
Sbias = [8*pi/180,8*pi/180,8*pi/180];        % Sun angle measurement bias
Sn    = [0.4*pi/180,0.4*pi/180,0.4*pi/180];  % Sun angle measurement noise
Torbit= 2*pi*sqrt((a0^3)/mu);                % Orbit Period
dt    = 0.5;                                 % Loop Period
[R0,V0] = COE2RV(e0,i0,Raan0,aop0,ta0,p0,mu);%(e,i,raan,aop,ta,p,mu)
s01   = [w0*pi/180  w0*pi/180 w0*pi/180 0.378 -0.378 0.756 0.378 R0(1,1) R0(2,1) R0(3,1) V0(1,1) V0(2,1) V0(3,1)];
t01   = 0:dt:30000;
[t,s] = ode23tb(@test27,t01,s01);
q     = [s(:,4) s(:,5) s(:,6) s(:,7)];
At    = quat2dcm(q);
[yaw, pitch, roll] = quat2angle(q);
SunIp = [cos(2*pi*(t01)/31557600);...
         cos(23.44*pi/180)*sin(2*pi*(t01)/31557600);...
         sin(23.44*pi/180)*sin(2*pi*(t01)/31557600)];
SunI  = SunIp';
gammaE1  = asin((RE)/(RE+alt));
for n= 1:numel(t01)
   alpha1(1,n) = acos(dot(s(n,8:10),-SunI(n,:)) / sqrt(s(n,8)^2+s(n,9)^2+s(n,10)^2));
if gammaE1  > alpha1(n)
    elipsoid(n) = 1;
else
    elipsoid(n) = 0;
end
phio(n,:)   = mod((real(dot(sign(-s(n,9)),acos((-s(n,8))/(sqrt(s(n,8)^2+s(n,9)^2)))))),2*pi);
tetao(n,:)  = acos((s(n,10))/(sqrt(s(n,8)^2+s(n,9)^2+s(n,10)^2)));
vs(n,:)     =  [cos(tetao(n,:)),0,-sin(tetao(n,:));0,1,0;sin(tetao(n,:)),0,cos(tetao(n,:))]...
               *[cos(phio(n,:)),sin(phio(n,:)),0;-sin(phio(n,:)),cos(phio(n,:)),0;0,0,1]...
               *[s(n,11);s(n,12);s(n,13)];  
sawo(n,:)   = mod((real(dot(sign(-vs(n,2)),acos((-vs(n,1))/(sqrt(vs(n,1)^2+vs(n,2)^2)))))),2*pi);
AOI         = [cos(sawo(n,:)),sin(sawo(n,:)),0;-sin(sawo(n,:)),cos(sawo(n,:)),0;0,0,1]*...
              [cos(tetao(n,:)),0,-sin(tetao(n,:));0,1,0;sin(tetao(n,:)),0,cos(tetao(n,:))]*...
              [cos(phio(n,:)),sin(phio(n,:)),0;-sin(phio(n,:)),cos(phio(n,:)),0;0,0,1];
AIO         = AOI';          
qd(n,:)     = (dcm2quat(AIO))';                              % q Desire Control algorithm
wd(n,:)     = (2*pi/Torbit)*AIO*[0;0;1];                     % w Desire Control algorithm
end
phr    = (2*pi*(t01)/86163.9)';
mres   = [0.005,0.005,0.005];
As     = 0.034;
rho    = 3*10^-9;
Cd     = 2.2;
Do     = 0.5*rho*Cd*As*((2*pi*a0/Torbit)^2)*[-1;0;0];
Di     = AIO*Do;
Db     = quat2dcm(q(1,:))*Di;
rcp    = [0.005;0.005;0.005];
Taerob = cross(rcp,Db);
for k=1:numel(t01)
lla(k,:)        = eci2lla([s(k,8) s(k,9) s(k,10)],datevec(datenum(time0)+(k-1)*dt*0.00001157405));
Bned(k,:)       = igrf((datenum(time0)+(k-1)*dt*0.00001157405),lla(k,1),lla(k,2),lla(k,3),'geodetic');
[Becef(k,1) Becef(k,2) Becef(k,3)]= ned2ecefv(Bned(k,1),Bned(k,2),Bned(k,3),lla(k,1),lla(k,2));
TEI             = [cos(phr(k,1)),sin(phr(k,1)),0;-sin(phr(k,1)),cos(phr(k,1)),0;0,0,1];
B(k,:)          = TEI^-1*Becef(k,:)';
Bb(k,:)         = quat2dcm([s(k,4) s(k,5) s(k,6) s(k,7)])*B(k,:)';
Tres(k,:)       = cross(mres,Bb(k,:));
Sunb(k,:)       = quat2dcm([s(k,4) s(k,5) s(k,6) s(k,7)])*SunI(k,:)';
Bm(k,:)         = Bb(k,:) + Bn + Bbias;
ABO             = quat2dcm([s(k,4) s(k,5) s(k,6) s(k,7)])*AIO;
AOB             = ABO';
Tg(k,:)         = cross(3*(nO^2)*AOB(:,3),I*AOB(:,3));
end
Smangle      = Sbias + Sn;
% qsunerror  =  angle2quat(Smangle(3), Smangle(2), Smangle(1));
Smp          = angle2dcm(Smangle(3), Smangle(2), Smangle(1))*Sunb';
Sm           = Smp';   
%% TRIAD method
for p = 1:numel(t01)
    n2(:,p)    = SunI(p,:)'/norm(SunI(p,:)');
    n1(:,p)    = B(p,:)'/norm(B(p,:)');
    b2(:,p)    = Sm(p,:)'/norm(Sm(p,:)');
    b1(:,p)    = Bm(p,:)'/norm(Bm(p,:)');
    %BT matrix
    bt2(:,p)   = cross(b1(:,p),b2(:,p))/norm(cross(b1(:,p),b2(:,p)));
    bt3(:,p)   = cross(b1(:,p),bt2(:,p));
    bt(:,:,p)  = [b1(:,p) bt2(:,p) bt3(:,p)];
    %NT matrix
    nt2(:,p)   = cross(n1(:,p),n2(:,p))/norm(cross(n1(:,p),n2(:,p)));
    nt3(:,p)   = cross(n1(:,p),nt2(:,p));
    nt(p,:,:)  = [transpose(n1(:,p));transpose(nt2(:,p));transpose(nt3(:,p))];
    %Calculate Roll,Pitch,Yaw measurement
    [yawm(p), pitchm(p), rollm(p)] = dcm2angle([b1(1,p)*n1(1,p)+bt2(1,p)*nt2(1,p)+bt3(1,p)*nt3(1,p),...
                  b1(1,p)*n1(2,p)+bt2(1,p)*nt2(2,p)+bt3(1,p)*nt3(2,p),...
                  b1(1,p)*n1(3,p)+bt2(1,p)*nt2(3,p)+bt3(1,p)*nt3(3,p);...
                  b1(2,p)*n1(1,p)+bt2(2,p)*nt2(1,p)+bt3(2,p)*nt3(1,p),...
                  b1(2,p)*n1(2,p)+bt2(2,p)*nt2(2,p)+bt3(2,p)*nt3(2,p),...
                  b1(2,p)*n1(3,p)+bt2(2,p)*nt2(3,p)+bt3(2,p)*nt3(3,p);...
                  b1(3,p)*n1(1,p)+bt2(3,p)*nt2(1,p)+bt3(3,p)*nt3(1,p),...
                  b1(3,p)*n1(2,p)+bt2(3,p)*nt2(2,p)+bt3(3,p)*nt3(2,p),...
                  b1(3,p)*n1(3,p)+bt2(3,p)*nt2(2,p)+bt3(3,p)*nt3(3,p)]);
end
%% Attitude Estimation
qk(1,:) = s(1,4:7);
wk(1,:) = s(1,1:3);
kd      = -0.2;
kp      = -0.02;
x(1,:)  = [qk(1,:) wk(1,:)];
sigmaq  = 0.005;
sigmaw  = 0.00003;
sigmab  = 3.82;
sigmas  = 3;
Pk      = eye(7)*0.10000;
for z = 1:numel(t01) -1
% Control Algorithm
qe(:,z)    = [ qd(z,4), qd(z,3),-qd(z,2),-qd(z,1);...
              -qd(z,3), qd(z,4), qd(z,1),-qd(z,2);...
               qd(z,2),-qd(z,1), qd(z,4),-qd(z,3);...
               qd(z,1), qd(z,2), qd(z,3), qd(z,4)]...
               *[x(z,1);x(z,2);x(z,3);x(z,4)];
we(:,z)    = [x(z,5);x(z,6);x(z,7)] - wd(z,:)';

%   propagate
    Omega    = [   0   , x(z,7),-x(z,6),x(z,5);...
                -x(z,7),   0   , x(z,5),x(z,6);...
                 x(z,6),-x(z,5),   0   ,x(z,7);...
                -x(z,5),-x(z,6),-x(z,7),  0  ]; 
    saw      = expm(0.5*Omega*dt);
    zeta     = [ x(z,4),-x(z,3), x(z,2);...
                 x(z,3), x(z,4),-x(z,1);...
                -x(z,2), x(z,1), x(z,4);...
                -x(z,1),-x(z,2),-x(z,3)];
    theta    = [0,x(z,7)*((I(2,2)-I(3,3))/I(1,1)),x(z,6)*((I(2,2)-I(3,3))/I(1,1));...
                x(z,7)*((I(3,3)-I(1,1))/I(2,2)),0,x(z,5)*((I(3,3)-I(1,1))/I(2,2));...
                x(z,5)*((I(1,1)-I(2,2))/I(3,3)),x(z,7)*((I(1,1)-I(2,2))/I(3,3)),0];
    Fp       = [0.5*Omega , 0.5*zeta;...
                zeros(3,4),  theta ];
    say      = expm(Fp*dt);
    Q        = [(sigmaq^2)*eye(4),zeros(4,3);...
                zeros(3,4),(sigmaw^2)*eye(3)];
    Pk       = say * Pk *say' + Q;
    Tc(:,z)  = kd*I*we(:,z)  + kp*I*qe(1:3,z);
    Td(:,z)  = Taerob + Tres(z,:)' + Tg(z,:)';
    dx5      = (Td(1,z) + Tc(1,z) + ( - x(z,6) * x(z,7) * ( I(3,3) - I(2,2) )))/I(1,1);
    dx6      = (Td(2,z) + Tc(2,z) + ( - x(z,5) * x(z,7) * ( I(1,1) - I(3,3) )))/I(2,2);
    dx7      = (Td(3,z) + Tc(3,z) + ( - x(z,5) * x(z,6) * ( I(2,2) - I(1,1) )))/I(3,3);
    x(z+1,:) = [(saw*x(z,1:4)')',x(z,5:7)+[dx5,dx6,dx7]*dt];

    
    
    
    
%   update
                        cb(z+1,:) = 0.5*(Bm(z+1,:) + B(z+1,:));
                        bb(z+1,:) = 0.5*(Bm(z+1,:) - B(z+1,:));
                        Hb        = [     0    , cb(z+1,3),-cb(z+1,2), bb(z+1,1);...
                                     -cb(z+1,3),     0    , cb(z+1,1), bb(z+1,2);...
                                      cb(z+1,2),-cb(z+1,1),    0   , bb(z+1,3);...
                                     -bb(z+1,1),-bb(z+1,2),-bb(z+1,3),   0   ];
                        Hbbar     = [Hb,zeros(4,3)];
                        Rb        = (sigmab^2)*eye(4);
                        Kb        = (Pk*Hbbar')/(Hbbar*Pk*Hbbar'+Rb);
                        x(z+1,:)  = (eye(7) - Kb*Hbbar)*x(z+1,:)';
                        Pk        = (eye(7) - Kb*Hbbar)*Pk*(eye(7) - Kb*Hbbar)' + Kb*Rb*Kb';
if elipsoid==0                        
cs(z+1,:) = 0.5*(Sm(z+1,:) + SunI(z+1,:));
bs(z+1,:) = 0.5*(Sm(z+1,:) - SunI(z+1,:));
Hs        = [    0     , cs(z+1,3),-cs(z+1,2), bs(z+1,1);...
             -cs(z+1,3),    0     , cs(z+1,1), bs(z+1,2);...
              cs(z+1,2),-cs(z+1,1),     0    , bs(z+1,3);...
             -bs(z+1,1),-bs(z+1,2),-bs(z+1,3),   0   ];
Hsbar     = [Hs,zeros(4,3)];
Rs        = (sigmas^2)*eye(4);
Ks        = (Pk*Hsbar')/(Hsbar*Pk*Hsbar'+Rs);
x(z+1,:)  = (eye(7) - Ks*Hsbar)*x(z+1,:)';
Pk        = (eye(7) - Ks*Hsbar)*Pk*(eye(7) - Ks*Hsbar)' + Ks*Rs*Ks';
end
end
qestimate    = x(:,1:4);
westimate    = x(:,5:7);
[yawestimate, pitchestimate, rollestimate] = quat2angle(qestimate);
qhat         = quatmultiply(quatconj(q),qestimate);
qdelta       = 180 - 2*acos(abs(qhat(:,4)))*180/pi;
%% Figures 
% figure
% hold all;
% [earth_x,earth_y,earth_z] = sphere;
% mesh (earth_x*RE,earth_y*RE,earth_z*RE);
% hold on
% plot3(s(:,8),s(:,9),s(:,10));
% grid on
% xlabel('X(km)');
% ylabel('Y(km)');
% zlabel('Z(km)')
% hold on
figure
subplot(3,1,1);
plot(t01,roll*180/pi,'r',t01,rollm*180/pi,':',t01,rollestimate*180/pi,'k--');
hold on
xlabel('time(sec)');
ylabel('roll(deg)');
legend('Roll Real','Roll Measurement','Roll Estimation')
grid on
subplot(3,1,2);
plot(t01,pitch*180/pi,'r',t01,pitchm*180/pi,':',t01,pitchestimate*180/pi,'k--');
hold on
xlabel('time(sec)');
ylabel('pitch(deg)');
legend('Pitch Real','Pitch Measurement','Pitch Estimation')
grid on
subplot(3,1,3);
plot(t01,yaw*180/pi,'r',t01,yawm*180/pi,':',t01,yawestimate*180/pi,'k--');
hold on
xlabel('time(sec)');
ylabel('yaw(deg)');
legend('Yaw Real','Yaw Measurement','Yaw Estimation')
grid on
% figure
% subplot(3,1,1);
% plot(t01-datenum(time0),rollm*180/pi);
% hold on
% xlabel('time(sec)');
% ylabel('roll measure(deg)');
% grid on
% subplot(3,1,2);
% plot(t01-datenum(time0),pitchm*180/pi);
% hold on
% xlabel('time(sec)');
% ylabel('pitch measure(deg)');
% grid on
% subplot(3,1,3);
% plot(t01-datenum(time0),yawm*180/pi);
% hold on
% xlabel('time(sec)');
% ylabel('yaw measure(deg)');
% grid on
figure
subplot(3,1,1);
plot(t01,s(:,1)*180/pi,t01,westimate(:,1)*180/pi);
hold on
xlabel('time(sec)');
ylabel('w_x(deg/sec)');
legend('W_x Real','W_x Estimation')
grid on
subplot(3,1,2);
plot(t01,s(:,2)*180/pi,t01,westimate(:,2)*180/pi);
hold on
xlabel('time(sec)');
ylabel('w_y(deg/sec)');
legend('W_y Real','W_y Estimation')
grid on
subplot(3,1,3);
plot(t01,s(:,3)*180/pi,t01,westimate(:,3)*180/pi);
hold on
xlabel('time(sec)');
ylabel('w_z(deg/sec)');
legend('W_z Real','W_z Estimation')
grid on
figure
plot(t01,qdelta,t01,elipsoid==1);
grid on
legend('deltaEstimation','Eclipse')