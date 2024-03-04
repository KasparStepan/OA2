clear all
close all
clc

c=1.42;
x2=0.65;
h=0.08;

H=0;
T=273.15-0.0065*H+15;
p=101325; %Pa, mozno vyuzit funkci v matlabu dle vysky letu, pro toto je zadan manualne 
M = 1.6;
ro=1.225;
v=336*M;
kappa=1.41;

AoA=[-3,-1,0,1,2,3,4,5];
deg2rad(AoA)

F=zeros(length(AoA),3);

for i = 1:length(AoA)
    delta1(i)=atan(h/(c-x2))-deg2rad(AoA(i));
    delta2(i)=atan(h/x2)+delta1(i)+deg2rad(AoA(i));
    delta3(i)=deg2rad(AoA(i));

    %% Prvni usek

    [p1(i),M1(i)]=solveCompression(delta1(i),M,p,kappa);
    F(i,1)=(p1(i)*sqrt((c-x2)^2+h^2));

    %% Druhz usek

    [p2(i),M2(i)]=solveExpansion(delta2(i),M1(i),p1(i),kappa);
    F(i,2)=(p2(i)*sqrt((x2)^2+h^2));
    disp("EV")
    %% Treti usek
    if delta3(i)==0
        p3(i)=p;
        M3(i)=M;
        F(i,3)=(p3(i)*c);
    elseif delta3(i)<0
        [p3(i),M3(i)]=solveExpansion(abs(delta3(i)),M,p,kappa);
        F(3)=(p3(i)*c);
        disp("EV")
    else
        [p3(i),M3(i)]=solveCompression(abs(delta3(i)),M,p,kappa);
        F(3)=(p3(i)*c);
    end


    %% Silová rovnováha
    b1=atan(h/(c-x2));
    b2=atan(h/x2);

    F1t(i) = p1(i)*(h^2+(c-x2)^2)^(1/2)*sin(b1); %t...rovnobezne s tetivou profilu 
    F1n(i) = -p1(i)*(h^2+(c-x2)^2)^(1/2)*cos(b1); %n...kolme na tetivu profilu 
    F2t(i) = -p2(i)*(h^2+x2^2)^(1/2)*sin(b2);
    F2n(i) = -p2(i)*(h^2+x2^2)^(1/2)*cos(b2);
    F3n(i) = p3(i)*c;
    
    N(i) = F1n(i)+F2n(i)+F3n(i); %sila kolmo na tetivu
    T(i) = F1t(i)+F2t(i); %sila rovnobezne s tetivou
    
    D(i) = N(i)*sin(deg2rad(AoA(i)))+T(i)*cos(deg2rad(AoA(i)));
    L(i) = N(i)*cos(deg2rad(AoA(i)))-T(i)*sin(deg2rad(AoA(i)));
    

    Mom(i)=F1t(i)*h/2+F2t(i)*h/2+F1n(i)*x2/2+F2n(i)*(x2/2-c/2);
    cL(i)=solveCoefficient(L(i),ro,v,c,"force");
    cD(i)=solveCoefficient(D(i),ro,v,c,"force");
    cM(i)=solveCoefficient(Mom(i),ro,v,c,"moment");

    fprintf("\nVýsledky na základě nelinearizované teorie\n")
    fprintf("Úhel náběhu: %d \n",AoA(i))
    fprintf("Vztlak profilu: %d N\n",L(i))
    fprintf("Odpor profilu: %d N\n",D(i))
    fprintf("Kroutící moment profilu: %d N/m\n",Mom(i))
    fprintf("Koeficient vztlaku profilu: %d \n",cL(i))
    fprintf("Koeficient odporu profilu: %d \n",cD(i))
    fprintf("Koeficient momentu profilu: %d \n",cM(i))
end

    x1 = 0;
    x0 = c/2;
    x2 = c-x2;
    Dx1 = c-x2;
    Dx2 = x2;
    Dyh1 = h;
    Dyh2 = -h;
    

for j=1:length(AoA)

    cll(j) = 4*deg2rad(AoA(j))/((M^2-1)^(1/2)); 
    cdl(j) = 4*deg2rad(AoA(j))^2/((M^2-1)^(1/2))+2/((M^2-1)^(1/2))*((Dyh1/Dx1)^2*Dx1/c+((Dyh2/Dx2)^2*Dx2/c));
    cml(j) = -4*deg2rad(AoA(j))/((M^2-1)^(1/2))*(1/2-x0/c)+2/((M^2-1)^(1/2))*((Dyh1/Dx1*(x1-x0)/c*Dx1/c)+(Dyh2/Dx2)*(x2-x0)/c*Dx2/c);
    
    fprintf("\nVýsledky na základě linearizované teorie\n")
    fprintf("Úhel náběhu: %d \n",AoA(j))
    fprintf("Koeficient vztlaku profilu: %d \n",cll(j))
    fprintf("Koeficient odporu profilu: %d \n",cdl(j))
    fprintf("Koeficient momentu profilu: %d \n",cml(j))
end

plot(AoA,cL)
hold on
plot(AoA,cll)
hold on
plot(AoA,cD)
hold on
plot(AoA,cdl)
hold on



function coefficient = solveCoefficient(force,ro,v,dimension,type)
if type == "force"
    coefficient=2*force/(ro*v^2*dimension);
elseif type == "moment"
    coefficient=2*force/(ro*v^2*dimension^2);
end
end

function [pressure, mach] = solveExpansion(d,M,p,kappa)

a=1;
b=3;
delta_i=1
error=10^-6;
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