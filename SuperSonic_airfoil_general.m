clear all
clc
close all

%% Okrajove podminky
pressure_atm=101325;
ro = 1.225;
T = 273.15;
Mach_number = 1.6;
kappa = 1.41;

error = 10^(-6)

%% Uhel nabehu
AoA = -5:5
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



top_airfoil_base = [0 top_surface_X 1; 0 top_surface_Y 0]
bottom_airfoil_base = [0 bottom_surface_X 1; 0 bottom_surface_Y 0]

if geometry_Validation(top_airfoil_base,bottom_airfoil_base) == true

    p = zeros(2,max([length(top_airfoil_base) length(bottom_airfoil_base)]),length(AoA));
    M = zeros(2,max([length(top_airfoil_base) length(bottom_airfoil_base)]),length(AoA));

    p(:,1,:) = pressure_atm;
    M(:,1,:) = Mach_number;

    L = zeros(1,length(AoA));
    D = zeros(1,length(AoA));
    Mom = zeros(1,length(AoA));

    Cl = zeros(1,length(AoA));
    Cd = zeros(1,length(AoA));
    Cm = zeros(1,length(AoA));

    for j = 1:length(AoA)

        top_airfoil = transformAirfoil(top_airfoil_base,AoA(j),chord);
        bottom_airfoil = transformAirfoil(bottom_airfoil_base, AoA(j),chord);

        [delta,normal] = getAngles(AoA(j),top_airfoil,bottom_airfoil)
        
        

        for i = 1:length(top_airfoil(1,:))-1
            if delta(1,i)>0
                [p(1,i+1,j), M(1,i+1,j)] = solveCompression(delta(1,i),M(1,i,j),p(1,i,j),kappa);
                
            elseif delta(1,i)<0
                
                [p(1,i+1,j), M(1,i+1,j)] = solveExpansion(abs(delta(1,i)),M(1,i,j),p(1,i,j),kappa);
            else
                p(1,i+1,j)=p(1,i,j);
                M(1,i+1,j)=M(1,i,j);
                

            end
        end
        for i = 1:length(bottom_airfoil(1,:))-1

            if delta(2,i)<0
                [p(2,i+1,j), M(2,i+1,j)] = solveCompression(abs(delta(2,i)),M(2,i,j),p(2,i,j),kappa);
                
            elseif delta(2,i)>0
                [p(2,i+1,j), M(2,i+1,j)] = solveExpansion(abs(delta(2,i)),M(2,i,j),p(2,i,j),kappa);
                
            else
                p(2,i+1,j)=p(2,i,j);
                M(2,i+1,j)=M(2,i,j);
                
            end

        end

    [L(j),D(j),Mom(j)] = solveForceMoment(p(:,:,j),normal, top_airfoil, bottom_airfoil);
    
    Cl(j) = solveCoefficient(L(j),ro,Mach_number,chord,"force");
    Cd(j) = solveCoefficient(D(j),ro,Mach_number,chord,"force");
    Cm(j) = solveCoefficient(Mom(j),ro,Mach_number,chord,"moment");

    fprintf("\nVýsledky na základě nelinearizované teorie\n")
    fprintf("Úhel náběhu: %d \n",AoA(j))
    fprintf("Vztlak profilu: %d N\n",L(j))
    fprintf("Odpor profilu: %d N\n",D(j))
    fprintf("Kroutící moment profilu: %d N/m\n",Mom(j))
    fprintf("Koeficient vztlaku profilu: %d \n",Cl(j))
    fprintf("Koeficient odporu profilu: %d \n",Cd(j))
    fprintf("Koeficient momentu profilu: %d \n",Cm(j))

    
    end
    [cll,cdl] = solveLinearTheory(AoA,Mach_number,top_airfoil_base,bottom_airfoil_base)
end

subplot(3,1,1)
plot(AoA,Cl)
grid on
title("Vztlaková čára")

subplot(3,1,2)
plot(AoA,Cd)
grid on
title("Odporová čára")

subplot(3,1,3)
plot(Cd,Cl)
grid on
title("Aerodynamická polára")




    


%plot(top_airfoil(1,:),top_airfoil(2,:))
%hold on
%plot(bottom_airfoil(1,:),bottom_airfoil(2,:))
%axis equal;
%grid on;







%% Definice funkci

function [changedSA] = transformAirfoil(SA,AoA,chord)

changedSA(1,:) = chord*(SA(2,:)*sin(deg2rad(AoA))+SA(1,:)*cos(deg2rad(AoA)));
changedSA(2,:) = chord*(SA(2,:)*cos(deg2rad(AoA))-SA(1,:)*sin(deg2rad(AoA)));


end

% Tato funnkce ověřuje jestli byla geometrie zadána správně, či nikoli
function [validation] = geometry_Validation(TA,BA)
    % Tato podmínka určuje zda počet X a Y souřadnic bodů je stejný pro
    % definici vrchní křivky profilu respektive spodní
    if (length(TA(1,:)) == length(TA(2,:))) & (length(BA(1,:)) == length(BA(2,:)))
        disp("Geometrie profilu je zadaná správně")
        validation = true;
    else
        disp("Geometrie profilu je zadaná špatně, nutno zkontrolovat správné zadání bodů")
        validation = false;
    end
end


% Tato funkce vypočítává z geometrie relativní úhly mezi jednotlivými
% hranami. Tyto úhly jsou potřebné pro výpočet kompresní nebo rázové vlny.
% Zároveň funkce vypočte normály křivek, které poslouží pro výpoet vztlaku
% a odporu z generovaných tlaků na hranách
function [delta, normal] = getAngles(AoA, TA, BA)

    for i = 1:(length(TA)-1)
        if i>1
            % Výpočet úhlu relativního vůči předešlé hrany, je to ypočítáno
            % s pomocí směrnice dané křivky, od které je odečten relativní
            % úhel předchozí křívky. 
            delta(1,i) =  atan((TA(2,i+1)-TA(2,i))/(TA(1,i+1)-TA(1,i))) -delta(1,i-1);
            
        else
            % Tato funkce vypočte relativní úhel pro první hranu bez
            % uvažování předešlé hrany.
            delta(1,i) =  atan((TA(2,i+1)-TA(2,i))/(TA(1,i+1)-TA(1,i)));
        end
        normal(1,i) = atan((TA(2,i+1)-TA(2,i))/(TA(1,i+1)-TA(1,i)));

    end
    % Viz komentáře v předchozím bloku for, v tomto případě je výpočet úhlů pro spodní stranu profil 
    for j = 1:(length(BA)-1)
        if j>1
           
            delta(2,j) = atan((BA(2,j+1)-BA(2,j))/(BA(1,j+1)-BA(1,j))-delta(2,j));
            
        else
            delta(2,j) = atan((BA(2,j+1)-BA(2,j))/(BA(1,j+1)-BA(1,j)));
           
        end
        normal(2,j) = atan((BA(2,j+1)-BA(2,j))/(BA(1,j+1)-BA(1,j)));

    end
end

% Funkce, která z řeší expanzní vlnu na základě relativního úhlu dané
% hrany, tlaku a machova čísla na předešlé hraně. Je využito metody
% bisekce.
function [pressure, mach] = solveExpansion(d,M,p,kappa)

a=1;
b=3;
delta_i=1;
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

% Funkce, která z řeší kompresní vlnu na základě relativního úhlu dané
% hrany, tlaku a machova čísla na předešlé hraně. Je využito metody
% bisekce.
function [pressure, mach] = solveCompression(d,M,p,kappa)
disp("RV");
a=deg2rad(30);
b=deg2rad(70);
delta_i=1;
error=10^-6;
while abs(delta_i)>error

    sigma_help=(a+b)/2;
    rad2deg(sigma_help);
    fun_sigma=(2*(M^2*sin(sigma_help)^2-1)/tan(sigma_help))/(2+M^2*(kappa+cos(2*sigma_help)));
    delta_i=(fun_sigma-d)/d;
    if delta_i<0
        a=sigma_help;
        rad2deg(sigma_help);
    else
        b=sigma_help;
        rad2deg(sigma_help);
    end
  
       
end
a;
b;
pressure = p*(2*kappa*M^2*sin(sigma_help)^2-(kappa-1))/(kappa+1);
mach = (((kappa-1)*M^2*sin(sigma_help)^2+2)/(2*kappa*M^2*(sin(sigma_help))^2-(kappa-1))/(sin(sigma_help-d))^2)^(1/2);

end

function [L,D,M] = solveForceMoment(p,direction, coordinates_top, coordinates_bottom)

L_top=[];
D_top=[];
L_bottom=[];
D_bottom=[];
Mom = 0;
for i = 1: length(coordinates_top(1,:))-1
    edge_length = sqrt((coordinates_top(1,i+1)-coordinates_top(1,i))^2+(coordinates_top(2,i+1)-coordinates_top(2,i))^2);
    L_top(i)= -1*p(1,i+1)*edge_length*cos(direction(1,i));
    D_top(i) = p(1,i+1)*edge_length*sin(direction(1,i));
    Mom = Mom - L_top(i)*((coordinates_top(1,i+1)+coordinates_top(1,i))/2-(coordinates_top(1,end)+coordinates_top(1,1))/2) + D_top(i)*((coordinates_top(2,i+1)+coordinates_top(2,i))/2-(coordinates_top(2,end)+coordinates_top(2,1))/2);
end    
for j = 1: length(coordinates_bottom(1,:))-1
    edge_length = sqrt((coordinates_bottom(1,j+1)-coordinates_bottom(1,j))^2+(coordinates_bottom(2,j+1)-coordinates_bottom(2,j))^2);
    L_bottom(j)= p(2,j+1)*edge_length*cos(direction(2,j));
    D_bottom(j) = -1*p(2,j+1)*edge_length*sin(direction(2,j));
    Mom = Mom - L_bottom(j)*((coordinates_bottom(1,j+1)+coordinates_bottom(1,j))/2-(coordinates_bottom(1,end)+coordinates_bottom(1,1))/2) + D_bottom(j)*((coordinates_bottom(2,j+1)+coordinates_bottom(2,j))/2-(coordinates_bottom(2,end)-coordinates_bottom(2,1))/2);
end
M = Mom;
L = sum(L_top)+sum(L_bottom);
D = sum(D_top)+sum(D_bottom);

end


function coefficient = solveCoefficient(force,ro,Mach,dimension,type)

v = 336*Mach;

if type == "force"
    coefficient=2*force/(ro*v^2*dimension);
elseif type == "moment"
    coefficient=2*force/(ro*v^2*dimension^2);
end
end


function [cll,cdl] = solveLinearTheory(AoA, M, TA, BA)

for i=1:length(AoA)

cll(i) = 4*deg2rad(AoA(i))/((M^2-1)^(1/2));
cdl(i) = 4*deg2rad(AoA(i))^2/((M^2-1)^(1/2));


for j=1:length(TA(1,:))-1

    cdl(i) = cdl(i)+2/((M^2-1)^(1/2))*(((TA(2,j+1)-TA(2,j))/(TA(1,j+1)-TA(1,j)))^2*(TA(2,j+1)-TA(2,j))/(TA(1,1)-TA(1,end)));
   

end

for j=1:length(BA(1,:))-1

    cdl(i) = cdl(i)+2/((M^2-1)^(1/2))*(((BA(2,j+1)-BA(2,j))/(BA(1,j+1)-BA(1,j)))^2*(BA(2,j+1)-BA(2,j))/(BA(1,1)-BA(1,end)));

end
end

end