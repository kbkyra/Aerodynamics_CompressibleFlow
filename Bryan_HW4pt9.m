clc
clear
close all

%% Intitial Comments

%{
Kyra Bryan
AEM 413-001, Spring 2021
Project 1: Problem 4.9

Problem 4.9:
Consider the intersection of two shocks of opposite families, as sketched in
Fig. 4.23. For M1 = 3, pl = 1 atm, theta2 = 20deg, and theta3 = 15deg, calculate the
pressure in regions 4 and 4', and the flow direction phi, behind the refracted
shocks.
%}

%% Inputs and Constants
M1 = 3;
p1 = 1; %atm
theta2 = 20; %deg
theta3 = 15; %deg

gamma = 1.4;

%% Calculations for Region 2 

%eqn 4.18
a = (1 + ((gamma-1)/2)*M1*M1)*tand(theta2); 
b = (1 + ((gamma+1)/2)*M1*M1)*tand(theta2);
temp = roots([a -(M1*M1 - 1) b 1]);
beta2_roots = atand(temp);
beta2 = min(beta2_roots(beta2_roots>0)); %find weak beta of the 3 roots (negative, strong shock, weak shock)

Mn1_2 = M1*sind(beta2); %eqn 4.7
p2p1 = 1 + ((2*gamma)/(gamma+1))*(Mn1_2*Mn1_2-1); %eqn 4.9

Mn2 = sqrt((Mn1_2*Mn1_2 + (2/(gamma-1)))/(((2*gamma)/(gamma-1))*Mn1_2*Mn1_2 - 1)); %eqn 4.10
M2 = Mn2/(sind(beta2-theta2)); %eqn 4.12

%% Calculations for Region 3

%eqn 4.18
a = (1 + ((gamma-1)/2)*M1*M1)*tand(theta3);
b = (1 + ((gamma+1)/2)*M1*M1)*tand(theta3);
temp = roots([a -(M1*M1 - 1) b 1]);
beta3_roots = atand(temp);
beta3 = min(beta3_roots(beta3_roots>0)); %find weak beta of the 3 roots (negative, strong shock, weak shock)

Mn1_3 = M1*sind(beta3); %eqn 4.7
p3p1 = 1 + ((2*gamma)/(gamma+1))*(Mn1_3*Mn1_3-1); %eqn 4.9

Mn3 = sqrt((Mn1_3*Mn1_3 + (2/(gamma-1)))/(((2*gamma)/(gamma-1))*Mn1_3*Mn1_3 - 1)); %eqn 4.10
M3 = Mn3/(sind(beta3-theta3)); %eqn 4.12

%% Calculations for Regions 4 and 4' and phi

theta4 = linspace(0,theta2,1000); %deg, this will be less than theta2
theta4prime = linspace(0,theta2,1000); %deg

for i = 1:length(theta4)
    
    %region 4 (from region 3)
    %eqn 4.18
    a = (1 + ((gamma-1)/2)*M3*M3)*tand(theta4(i));
    b = (1 + ((gamma+1)/2)*M3*M3)*tand(theta4(i));
    temp = roots([a -(M3*M3 - 1) b 1]);
    beta4_roots = atand(temp);
    beta4 = min(beta4_roots(beta4_roots>0)); %find weak beta of the 3 roots (negative, strong shock, weak shock)

    Mn4 = M3*sind(beta4); %eqn 4.7
    p4p3 = 1 + ((2*gamma)/(gamma+1))*(Mn4*Mn4-1); %eqn 4.9
    p4 = p4p3*p3p1*p1;

    %region 4' (from region 2)
    for j = 1:length(theta4prime)
        %eqn 4.18
        a = (1 + ((gamma-1)/2)*M2*M2)*tand(theta4prime(j));
        b = (1 + ((gamma+1)/2)*M2*M2)*tand(theta4prime(j));
        temp = roots([a -(M2*M2 - 1) b 1]);
        beta4prime_roots = atand(temp);
        beta4prime = min(beta4prime_roots(beta4prime_roots>0)); %find weak beta of the 3 roots (negative, strong shock, weak shock)

        Mn4prime = M2*sind(beta4prime); %eqn 4.7
        p4primep2 = 1 + ((2*gamma)/(gamma+1))*(Mn4prime*Mn4prime-1); %eqn 4.9
        p4prime = p4primep2*p2p1*p1;

        %compare to see if p4 ~= p4prime
        if (p4-p4prime) > 0.0001
            phi = theta4(i)-theta3;
            p4final = p4;
            p4primefinal = p4prime;
        end 
    end
end 

%% Output 

fprintf('Pressure in region 4 = %.4f atm\n', p4final);
fprintf('Pressure in region 4'' = %.4f atm\n', p4primefinal);
fprintf('Flow direction phi = %.4f deg\n', phi); 

fprintf('\nKyra Bryan''s AEM413 Project: 4.9 script complete. --------------------------------------------------\n\n');