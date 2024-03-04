%% Okrajove podminky
p=101325;
ro = 101325;
T = 288.15;
M = 1.5;
kappa = 1.4;

error = 10^(-6)

%% Uhel nabehu
AoA = [5]
chord = 1.42;

%% Geometrie profilu
% Geometrie profilu je normovana na 1 metr delky, je proto nutne zadavat
% polohy bodů v poměru k délce tětivy, která má délku 1m

%Náběžná hrana a odtoková hrana jsou definovány jako [0,0] a [1,0]

top_surface_X = [0.65];
top_surface_Y = [0.3];

bottom_surface_X = [];
bottom_surface_Y = [];


%% Výpočet
%Vytvoření finální geometrie



top_airfoil = [0 top_surface_X 1; 0 top_surface_Y 0]
bottom_airfoil = [0 bottom_surface_X 1; 0 bottom_surface_Y 0]

if geometry_Validation(top_airfoil,bottom_airfoil) == true

top_airfoil = transformAirfoil(top_airfoil,AoA,chord);
bottom_airfoil = transformAirfoil(bottom_airfoil, AoA,chord);



delta = getAngles(AoA,top_airfoil,bottom_airfoil)


for i = 2:length(top_airfoil(1,:))
    if delta(1,i)>0
        [p(i), M(i)] = solveCompression(delta(1,i),M,p,kappa)
    elseif delta(1,i)<0
        [p(i), M(i)] = solveExpansion(delta(1,i),M,p,kappa)
    else
        p(i)=p(i-1)
        M(i)=M(i-1)



end
end
end







%% Definice funkci

function [changedSA] = transformAirfoil(SA,AoA,chord)

changedSA(1,:) = chord*(SA(2,:)*sin(deg2rad(AoA))+SA(1,:)*cos(deg2rad(AoA)));
changedSA(2,:) = chord*(SA(2,:)*cos(deg2rad(AoA))-SA(1,:)*sin(deg2rad(AoA)));


end


function [validation] = geometry_Validation(TA,BA)
    if (length(TA(1,:)) == length(TA(2,:))) & (length(BA(1,:)) == length(BA(2,:)))
        disp("Geometrie profilu je zadaná správně")
        validation = true;
    else
        disp("Geometrie profilu je zadaná špatně, nutno zkontrolovat správné zadání bodů")
        validation = false;
    end
end

function [direction] = getAngles(AoA, TA, BA)
disp(TA)
disp(BA)
    for i = 1:(length(TA)-1)
        if i>1
            direction(1,i) =  rad2deg(atan((TA(2,i+1)-TA(2,i))/(TA(1,i+1)-TA(1,i)))) -direction(1,i-1)
            disp("dir2")
            
        else
            direction(1,i) =  rad2deg(atan((TA(2,i+1)-TA(2,i))/(TA(1,i+1)-TA(1,i))))
            disp("dir1")
        end

    end

    for j = 1:(length(BA)-1)
        if j>1
            direction(2,j) = rad2deg(atan((BA(2,j+1)-BA(2,j))/(BA(1,j+1)-BA(1,j)))-direction(2,j))
        else
            direction(2,j) = rad2deg(atan((BA(2,j+1)-BA(2,j))/(BA(1,j+1)-BA(1,j))))
        end

    end
end

function [pressure, mach] = solveExpansion(d,M,p,kappa)

a=1;
b=3;
delta_i=1
error=10^-3;
v1=((kappa+1)/(kappa-1))^(1/2)*atan(((kappa-1)/(kappa+1)*(M^2-1))^(1/2))-atan((M^2-1)^(1/2));
v2=v1+d;
while abs(delta_i)>error
    mach_help=(a+b)/2;
    fun_mach=((kappa+1)/(kappa-1))^(1/2)*atan(((kappa-1)/(kappa+1)*(mach_help^2-1))^(1/2))-atan((mach_help^2-1)^(1/2));
    delta_i=(fun_mach-v2)/v2;
    
    if delta_i<0
        a=mach_help;
    else
        b=mach_help;
    end

end
a;
b;
mach = mach_help;
pressure = ((1+(kappa-1)/2*M^2)/(1+(kappa-1)/2*mach^2))^(kappa/(kappa-1))*p;

end

function [pressure, mach] = solveCompression(d,M,p,kappa)
a=deg2rad(30);
b=deg2rad(70);
delta_i=1;
error=10^-3;
while abs(delta_i)>error

    sigma_help=(a+b)/2
    rad2deg(sigma_help)
    fun_sigma=(2*(M^2*sin(sigma_help)^2-1)/tan(sigma_help))/(2+M^2*(kappa+cos(2*sigma_help)));
    delta_i=(fun_sigma-d)/d;
    if delta_i<0
        a=sigma_help;
    else
        b=sigma_help;
    end
  
       
end
a;
b;
pressure = p*(2*kappa*M^2*sin(sigma_help)^2-(kappa-1))/(kappa+1);
mach = (((kappa-1)*M^2*sin(sigma_help)^2+2)/(2*kappa*M^2*(sin(sigma_help))^2-(kappa-1))/(sin(sigma_help-d))^2)^(1/2);

end

function [L,D] = solveForce(p,direction, Pstart, Pend, surface)
length = sqrt((Pend(1)-Pstart(1))^2+(Pend(2)-Pstart(2))^2);

if surface == "top"
    L= p*length*cos(direction)
    D = p*length*sin(dirextio)
elseif surface == "bottom"
    L= -1*p*length*cos(direction)
    D = p*length*sin(dirextio)

end

end

function [M] = solveMoment(p,direction, Pstart, Pend, surface,Pref)

[L, D] = solveForce(p,direction, Pstart, Pend, surface)

Ml = L*((Pstart(1)+Pend(1))/2 - Pref(1))
Md = D*((Pstart(1)+Pend(1))/2 - Pref(2))

M = Ml + Md

end

function coefficient = solveCoefficient(force,ro,v,dimension,type)
if type == "force"
    coefficient=2*force/(ro*v^2*dimension);
elseif type == "moment"
    coefficient=2*force/(ro*v^2*dimension^2);
end
end
