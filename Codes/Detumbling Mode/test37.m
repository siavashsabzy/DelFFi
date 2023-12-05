function ds = test37(t,s)
I=[0.017,0,0;0,0.055,0;0,0,0.055];
time0=[2014 09 14 12 00 00];   % Satellite Start Mission
mu=398576.0576;
RE=6400;
J2=1.08263*10^-3;
alt=200;                         %  altitude in km
a0       = (alt+RE);             %  a    = semi-major axis, km
n        = sqrt(mu/(a0^3));      % Orbit angular Rate
Bbias = [500;500;500]*10^-9;     %  magnetic field measurement bias
Bn    = [150;150;150]*10^-9;     %  magnetic field measurement noise
phr   = (2*pi*(t)/86163.9)';
ds    = zeros(13,1);
lla   = eci2lla([s(8) s(9) s(10)],datevec(datenum(time0)+t*0.00001157405));
Bned  = igrf((datenum(time0)+t*0.00001157405),lla(1),lla(2),lla(3),'geod');
[Becef(1) Becef(2) Becef(3)]= ned2ecefv(Bned(1),Bned(2),Bned(3),lla(1),lla(2));
TEI             = [cos(phr(1)),sin(phr(1)),0;-sin(phr(1)),cos(phr(1)),0;0,0,1];
B               = TEI^-1*Becef';
A     = quat2dcm([s(4) s(5) s(6) s(7)]); %Inertial to body Transformation
phio  = mod((real(dot(sign(-s(9)),acos((-s(8))/(sqrt(s(8)^2+s(9)^2)))))),2*pi);
tetao = acos((s(10))/(sqrt(s(8)^2+s(9)^2+s(10)^2)));
vs    =  [cos(tetao),0,-sin(tetao);0,1,0;sin(tetao),0,cos(tetao)]...
        *[cos(phio),sin(phio),0;-sin(phio),cos(phio),0;0,0,1]*[s(11);s(12);s(13)];  
sawo  = mod((real(dot(sign(-vs(2)),acos((-vs(1))/(sqrt(vs(1)^2+vs(2)^2)))))),2*pi);
AIO   = [cos(sawo),sin(sawo),0;-sin(sawo),cos(sawo),0;0,0,1]*...
        [cos(tetao),0,-sin(tetao);0,1,0;sin(tetao),0,cos(tetao)]*...
        [cos(phio),sin(phio),0;-sin(phio),cos(phio),0;0,0,1];
As    = 0.034;
rho   = 3*10^-9;
Cd    = 2.2;
Torbit=2*pi*sqrt((a0^3)/mu);
Do    = 0.5*rho*Cd*As*((2*pi*a0/Torbit)^2)*[-1;0;0];
Di    = AIO*Do;
Db    = A*Di;
rcp   = [0.005;0.005;0.005];
Taerob = cross(rcp,Db);
Bb    =  A * B;
Bm    =  Bb + Bn + Bbias;
Bdot  =  cross(Bm,[s(1);s(2);s(3)]);
mres  = [0.005;0.005;0.005];
Tres  = cross(mres,Bb);
ABO   = quat2dcm([s(4) s(5) s(6) s(7)])*AIO;
AOB   = ABO';
kb    = 28000;
Tb    = cross((-kb*Bdot)/5,(Bm));               %  Case E ([t]actuate/[t]total = 1/5) so [Tb]eff = Tb/5  
Tg    =  cross(3*(n^2)*AOB(:,3),I*AOB(:,3));    %  GravityGradientTorque due to Moment of inertia
Tc1   =  Tb(1); Td1   = Taerob(1,1) + Tres(1,1) + Tg(1,1) ;
Tc2   =  Tb(2); Td2   = Taerob(2,1) + Tres(2,1) + Tg(2,1);
Tc3   =  Tb(3); Td3   = Taerob(3,1) + Tres(3,1) + Tg(3,1) ;
ds(1) = (Tc1 + Td1 + ( - s(2) * s(3) * ( I(3,3) - I(2,2) )))/I(1,1);
ds(2) = (Tc2 + Td2 + ( - s(1) * s(3) * ( I(1,1) - I(3,3) )))/I(2,2);
ds(3) = (Tc3 + Td3 + ( - s(1) * s(2) * ( I(2,2) - I(1,1) )))/I(3,3);
ds(4) = 0.5*( s(3)*s(5) - s(2)*s(6) + s(1)*s(7));
ds(5) = 0.5*(-s(3)*s(4) + s(1)*s(6) + s(2)*s(7));
ds(6) = 0.5*( s(2)*s(4) - s(1)*s(5) + s(3)*s(7));
ds(7) = 0.5*(-s(1)*s(4) - s(2)*s(5) - s(3)*s(6));
ds(8) = s(11);
ds(9) = s(12);
ds(10)= s(13);
ds(11)= (-mu/(sqrt(s(8)^2+s(9)^2+s(10)^2)^3))*s(8) + (3*mu*J2*(RE^2)/(2*(sqrt(s(8)^2+s(9)^2+s(10)^2)^5)))...
    *(5*((s(10)^2/(sqrt(s(8)^2+s(9)^2+s(10)^2)^2))-1))*s(8);
ds(12)= (-mu/(sqrt(s(8)^2+s(9)^2+s(10)^2)^3))*s(9) + (3*mu*J2*(RE^2)/(2*(sqrt(s(8)^2+s(9)^2+s(10)^2)^5)))...
    *(5*((s(10)^2/(sqrt(s(8)^2+s(9)^2+s(10)^2)^2))-1))*s(9);
ds(13)= (-mu/(sqrt(s(8)^2+s(9)^2+s(10)^2)^3))*s(10) + (3*mu*J2*(RE^2)/(2*(sqrt(s(8)^2+s(9)^2+s(10)^2)^5)))...
    *(5*((s(10)^2/(sqrt(s(8)^2+s(9)^2+s(10)^2)^2))-1)*s(10) - 2*s(10));