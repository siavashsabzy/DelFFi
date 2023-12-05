clear
clc
time0 = [2018 09 14 12 00 00];               %  Satellite Start Mission
I     = [0.017,0,0;0,0.055,0;0,0,0.055];     %  Satellite Moment of Inertia Matrix
mu    = 398576.0576;                         %  Gravitational Coe.
RE    = 6371;                                %  in km
alt   = 380;                                 %  altitude in km
a0    = (alt+RE);                            %  semi-major axis, km
e0    = 0;                                   %  eccentricity
i0    = 98.6;                                %  inclination, deg
Raan0 = -15;                                 %  Right ascension of the ascending node, deg
aop0  = 0;                                   %  argument of periapse, deg                           
ta0   = 0;                                   %  ta   = true anomaly, deg
p0    = a0*(1-e0^2);                         %  p0
nO    = sqrt(mu/(a0^3));                     %  Orbit angular Rate
w0    = 10;                                  %  degree per second
Bbias = [500,500,500]*10^-9;                 %  magnetic field measurement bias
Bn    = [150,150,150]*10^-9;                 %  magnetic field measurement noise
Torbit= 2*pi*sqrt((a0^3)/mu);                %  Orbit Period
dt    = 0.5;                                 %  Loop Period
[R0,V0] = COE2RV(e0,i0,Raan0,aop0,ta0,p0,mu);%(e,i,raan,aop,ta,p,mu)
s01=[w0*pi/180  w0*pi/180 w0*pi/180 0.378 -0.378 0.756 0.378 R0(1,1) R0(2,1) R0(3,1) V0(1,1) V0(2,1) V0(3,1)];
t01=0:dt:70;
[t,s]=ode23tb(@test37,t01,s01);
q = [s(:,4) s(:,5) s(:,6) s(:,7)];
At = quat2dcm(q);
[yaw, pitch, roll] = quat2angle(q);
phr = (2*pi*(t01)/86163.9)';
for k=1:numel(t01)
lla(k,:)        = eci2lla([s(k,8) s(k,9) s(k,10)],datevec(datenum(time0)+(k-1)*dt*0.00001157405));
Bned(k,:)       = igrf((datenum(time0)+(k-1)*dt*0.00001157405),lla(k,1),lla(k,2),lla(k,3),'geodetic');
[Becef(k,1) Becef(k,2) Becef(k,3)]= ned2ecefv(Bned(k,1),Bned(k,2),Bned(k,3),lla(k,1),lla(k,2));
TEI             = [cos(phr(k,1)),sin(phr(k,1)),0;-sin(phr(k,1)),cos(phr(k,1)),0;0,0,1];
B(k,:)          = TEI^-1*Becef(k,:)';
Bb(k,:)         = quat2dcm([s(k,4) s(k,5) s(k,6) s(k,7)])*B(k,:)';
Bm(k,:)         = Bb(k,:) + Bn + Bbias;
end
for z=2:numel(t01)
   Bdot(z,:)     = (Bb(z,:) - Bb(z-1,:))/dt; 
   Bdot(1,:)     = Bb(1,:)/dt; 
end

%%
figure
subplot(3,1,1);
plot(t01,s(:,1)*180/pi);
hold on
xlabel('time(sec)');
ylabel('w_x(deg/sec)');
grid on
subplot(3,1,2);
plot(t01,s(:,2)*180/pi);
hold on
xlabel('time(sec)');
ylabel('w_y(deg/sec)');
grid on
subplot(3,1,3);
plot(t01,s(:,3)*180/pi);
hold on
xlabel('time(sec)');
ylabel('w_z(deg/sec)');
grid on
figure
plot(t01,Bm(:,1),t01,Bm(:,2),t01,Bm(:,3))
figure
plot(t01,s(:,1)*180/pi,t01,s(:,2)*180/pi,t01,s(:,3)*180/pi)