%{
a = semi-major axis, m
e = eccentricity
i = inclination, deg
raan = big omega = right ascension of the ascending node, deg
aop = little omega = argument of periapse, deg
v = true anomaly, deg
R = [Ri,Rj,Rk] vector, m, ECI
V = [Vi,Vj,Vk] vector, m, ECI
mu = G*M(centralbody)
see algorithm 10. p. 118
*NOTE: capital variables are vectors; uncapitalized versions are their
norms
%}
function [R,V] = COE2RV(e,i,raan,aop,ta,p,mu)
    %0. constants and helper functions
    I = [1,0,0];
    J = [0,1,0];
    K = [0,0,1];
    syms rad deg
    r2d = @(rad) rad*180.0/pi;
    d2r = @(deg) deg*pi/180.0;
    
    %1. orbit type
    circle = 0;
    ellipse = 0;
    parabola = 0;
    hyperbola = 0;
    equatorial = 0;
    inclined = 0;
    if e == 0
        circle = 1;
    elseif e < 1
        ellipse = 1;
    elseif e == 1
        parabola = 1;
    elseif e > 1
        hyperbola = 1;
    end
    if i == 0 equatorial = 1; else inclined = 1; end;
    
    %2. fixes
    if circle && equatorial; aop = 0.0; raan = 0.0; end;
    if circle && inclined; aop = 0.0; end;
    if ellipse && equatorial; raan = 0.0; end;
    
    %3. R, V in PQW coords
    R_PQW1 = p*cosd(ta)/(1+e*cosd(ta));
    R_PQW2 = p*sind(ta)/(1+e*cosd(ta));
    R_PQW3 = 0;
    R_PQW = [R_PQW1; R_PQW2; R_PQW3];
    
    V_PQW1 = -1*sqrt(mu/p)*sind(ta);
    V_PQW2 = sqrt(mu/p)*(e+cosd(ta));
    V_PQW3 = 0;
    V_PQW = [V_PQW1; V_PQW2; V_PQW3];
    
    %3. transofmration matrix
    % column 1
    a11 = cosd(raan)*cosd(aop) - sind(raan)*sind(aop)*cosd(i);
    a21 = sind(raan)*cosd(aop) + cosd(raan)*sind(aop)*cosd(i);
    a31 = sind(aop)*sind(i);
    % column 2
    a12 = -1*cosd(raan)*sind(aop) - sind(raan)*cosd(aop)*cosd(i);
    a22 = -1*sind(raan)*sind(aop) + cosd(raan)*cosd(aop)*cosd(i);
    a32 = cosd(aop)*sind(i);
    % column 3
    a13 = sind(raan)*sind(i);
    a23 = -1*cos(raan)*sind(i);
    a33 = cosd(i);
    % final matrix
    T = [a11, a12, a13; a21, a22, a23; a31, a32, a33];
    
    %4. transform to R_IJK, V_IJK
    R = T*R_PQW;
    V = T*V_PQW;
end