

%% Okrajove podminky
press=101325;
ro = 1.225;
T = 273.15;
Mach = 1.6;
kappa = 1.41;

error = 10^(-6)

%% Uhel nabehu
AoA = [0]
chord = 1.42;

%% Geometrie profilu
% Geometrie profilu je normovana na 1 metr delky, je proto nutne zadavat
% polohy bodů v poměru k délce tětivy, která má délku 1m

%Náběžná hrana a odtoková hrana jsou definovány jako [0,0] a [1,0]

top_surface_X = [(1.42-0.65)/1.42];
top_surface_Y = [0.08/1.42];

bottom_surface_X = [];
bottom_surface_Y = [];


%% Výpočet
%Vytvoření finální geometrie



top_airfoil = [0 top_surface_X 1; 0 top_surface_Y 0]
bottom_airfoil = [0 bottom_surface_X 1; 0 bottom_surface_Y 0]

if geometry_Validation(top_airfoil,bottom_airfoil) == true

top_airfoil = transformAirfoil(top_airfoil,AoA,chord);
bottom_airfoil = transformAirfoil(bottom_airfoil, AoA,chord);



[delta,normal] = getAngles(AoA,top_airfoil,bottom_airfoil)
p(1,1) = press
p(2,1) = press
M(1,1) = Mach
M(2,1) = Mach

for i = 1:length(top_airfoil(1,:))-1
    if delta(1,i)>0
        [p(1,i+1), M(1,i+1)] = solveCompression(delta(1,i),M(1,i),p(1,i),kappa)
    elseif delta(1,i)<0
        [p(1,i+1), M(1,i+1)] = solveExpansion(abs(delta(1,i)),M(1,i),p(1,i),kappa)
    else
        p(1,i+1)=p(1,i)
        M(1,i+1)=M(1,i)

    end
for i = 1:length(bottom_airfoil(1,:))-1

    if delta(2,i)<0
        [p(2,i+1), M(2,i+1)] = solveCompression(abs(delta(2,i)),M(2,i),p(2,i),kappa)
    elseif delta(2,i)>0
        [p(2,i+1), M(2,i+1)] = solveExpansion(abs(delta(2,i)),M(2,i),p(2,i),kappa)
    else
        p(2,i+1)=p(2,i)
        M(2,i+1)=M(2,i)
    end



end
end
end



    


plot(top_airfoil(1,:),top_airfoil(2,:))
hold on
plot(bottom_airfoil(1,:),bottom_airfoil(2,:))
axis equal;
grid on;

[L,D] = solveForce(p,normal, top_airfoil, bottom_airfoil)





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

function [delta, normal] = getAngles(AoA, TA, BA)
disp(TA)
disp(BA)
    for i = 1:(length(TA)-1)
        if i>1
            delta(1,i) =  atan((TA(2,i+1)-TA(2,i))/(TA(1,i+1)-TA(1,i))) -delta(1,i-1)
            
        else
            delta(1,i) =  atan((TA(2,i+1)-TA(2,i))/(TA(1,i+1)-TA(1,i)))
        end
        normal(1,i) = atan((TA(2,i+1)-TA(2,i))/(TA(1,i+1)-TA(1,i)))

    end

    for j = 1:(length(BA)-1)
        if j>1
            delta(2,j) = atan((BA(2,j+1)-BA(2,j))/(BA(1,j+1)-BA(1,j))-delta(2,j))
            
        else
            delta(2,j) = atan((BA(2,j+1)-BA(2,j))/(BA(1,j+1)-BA(1,j)))
           
        end
        normal(2,j) = atan((BA(2,j+1)-BA(2,j))/(BA(1,j+1)-BA(1,j)))

    end
end

function [pressure, mach] = solveExpansion(d,M,p,kappa)

a=1;
b=3;
delta_i=1
error=10^-5;

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
disp("RV")
a=deg2rad(30);
b=deg2rad(70);
delta_i=1;
error=10^-6;
while abs(delta_i)>error

    sigma_help=(a+b)/2
    rad2deg(sigma_help)
    fun_sigma=(2*(M^2*sin(sigma_help)^2-1)/tan(sigma_help))/(2+M^2*(kappa+cos(2*sigma_help)));
    delta_i=(fun_sigma-d)/d;
    if delta_i<0
        a=sigma_help;
        rad2deg(sigma_help)
    else
        b=sigma_help;
        rad2deg(sigma_help)
    end
  
       
end
a;
b;
pressure = p*(2*kappa*M^2*sin(sigma_help)^2-(kappa-1))/(kappa+1);
mach = (((kappa-1)*M^2*sin(sigma_help)^2+2)/(2*kappa*M^2*(sin(sigma_help))^2-(kappa-1))/(sin(sigma_help-d))^2)^(1/2);

end

function [L,D] = solveForce(p,direction, coordinates_top, coordinates_bottom)

L_top=[]
D_top=[]
L_bottom=[]
D_bottom=[]
for i = 1: length(coordinates_top(1,:))-1
    edge_length = sqrt((coordinates_top(1,i+1)-coordinates_top(1,i))^2+(coordinates_top(2,i+1)-coordinates_top(2,i))^2);
    L_top(i)= -1*p(1,i+1)*edge_length*cos(direction(1,i))
    D_top(i) = p(1,i+1)*edge_length*sin(direction(1,i))
end    
for j = 1: length(coordinates_bottom(1,:))-1
    edge_length = sqrt((coordinates_bottom(1,j+1)-coordinates_bottom(1,j))^2+(coordinates_bottom(2,j+1)-coordinates_bottom(2,j))^2);
    L_bottom(j)= 1*p(2,j+1)*edge_length*cos(direction(2,j))
    D_bottom(j) = p(2,j+1)*edge_length*sin(direction(2,j))
end

L = sum(L_top)+sum(L_bottom);
D = sum(D_top)+sum(D_bottom);

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
