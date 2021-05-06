clc
clear
close all

%% Initial Comments

%{
Kyra Bryan
AEM 413-001, Spring 2021
Project 2: Problem 4.17

Problem 4.17:
Consider a flat plate with a chord length (from leading to trailing edge) of
1 m. The free-stream flow properties are M1 = 3, pl = 1 atm, and T1 =
270 K. Using shock-expansion theory, tabulate and plot on graph paper these
properties as functions of angle of attack from 0 to 30deg (use increments of 5deg):
a. Pressure on the top surface
b. Pressure on the bottom surface
c. Temperature on the top surface
d. Temperature on the bottom surface
e. Lift per unit span
f. Drag per unit span
g. Lift/drag ratio

Reference: Section 4.15 and corresponding examples
%}

%% Inputs, Constants

chord = 1; %m
M1 = 3;
P1 = 1; %atm
T1 = 270; %K

gamma = 1.4;

%% Calculations 

alpha = 0:5:30;
Ptop = zeros(length(alpha),1);
Ttop = zeros(length(alpha),1);
Pbot = zeros(length(alpha),1);
Tbot = zeros(length(alpha),1);
lift = zeros(length(alpha),1);
drag = zeros(length(alpha),1);
ld = zeros(length(alpha),1);

for i = 1:length(alpha)
      
    %% P, T - Top Surface
    
    nu1 = sqrt((gamma+1)/(gamma-1))*atand(sqrt(((gamma-1)*(M1^2-1))/(gamma+1))) - atand(sqrt(M1^2-1)); %eqn 4.44
    nu2 = alpha(i) + nu1;
    
    syms M
    eqn = nu2 == (sqrt((gamma+1)/(gamma-1))*atand(sqrt(((gamma-1)*(M^2-1))/(gamma+1)))) - atand(sqrt(M^2-1)); %eqn 4.44
    M2 = solve(eqn,M); 
    M2 = abs(M2); %correction for solve() returning the negative of the sqrt() sometimes
    
    Po1P1 = (1+((gamma-1)/2)*M1*M1)^(gamma/(gamma-1)); %eqn 3.30
    Po2P2 = (1+((gamma-1)/2)*M2*M2)^(gamma/(gamma-1)); %eqn 3.30
    Po2P1 = 1; %constant
    Ptop(i) = (1/Po2P2)*Po2P1*Po1P1*P1;
    
    To1T1 = 1 + ((gamma-1)/2)*M1*M1; %eqn 3.28
    To2T2 = 1 + ((gamma-1)/2)*M2*M2; %eqn 3.28
    To2T1 = 1; %constant
    Ttop(i) = (1/To2T2)*To2T1*To1T1*T1;
    
    %% P, T - Bottom Surface
    
    %eqn 4.18
    a = (1 + ((gamma-1)/2)*M1*M1)*tand(alpha(i));
    b = (1 + ((gamma+1)/2)*M1*M1)*tand(alpha(i));
    temp = roots([a -(M1*M1 - 1) b 1]);
    beta_roots = atand(temp);
    beta = min(beta_roots(beta_roots>0)); %find weak beta of the 3 roots (negative, strong shock, weak shock)

    Mn1 = M1*sind(beta); %eqn 4.7
    P3P1 = 1 + ((2*gamma)/(gamma+1))*(Mn1*Mn1-1); %eqn 4.9
    Pbot(i) = P3P1*P1;
    
    rho3rho1 = ((gamma+1)*Mn1*Mn1)/((gamma-1)*Mn1*Mn1+2); %eqn 4.8
    T3T1 = P3P1*(1/rho3rho1); %eqn 4.11
    Tbot(i) = T3T1*T1;
   
    %% L, D
    
    lift(i) = (Pbot(i)*101325-Ptop(i)*101325)*chord*cosd(alpha(i)); %converting P's atm to N/m^2 %ex 4.15 %N/m
    drag(i) = (Pbot(i)*101325-Ptop(i)*101325)*chord*sind(alpha(i)); %converting P's atm to N/m^2 %ex 4.15 %N/m
    ld(i) = lift(i)/drag(i);
    
end 

master_array = [alpha' Ptop Pbot Ttop Tbot lift drag ld];

%% Plotting, Output

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/20 scrsz(3)/20 scrsz(3)/1.25 scrsz(4)/1.25]) 

subplot(2,2,1)
plot(alpha,Ptop,'DisplayName','Pressure on top surface');
hold on
grid on
plot(alpha,Pbot,'DisplayName','Pressure on bottom surface');
xlabel("Angle of Attack (deg)");
ylabel("Pressure (atm)");
title("Fig 1a: Pressure vs. Angle of Attack");
legend('Location','best')

subplot(2,2,2)
plot(alpha,Ttop,'DisplayName','Temperature on top surface');
hold on
grid on
plot(alpha,Tbot,'DisplayName','Temperature on bottom surface');
xlabel("Angle of Attack (deg)");
ylabel("Temperature (K)");
title("Fig 1b: Temperature vs. Angle of Attack");
legend('Location','best')

subplot(2,2,3)
plot(alpha,lift,'DisplayName','Lift');
hold on
grid on
plot(alpha,drag,'DisplayName','Drag');
xlabel("Angle of Attack (deg)");
ylabel("Force Per Unit Span (N/m)");
title("Fig 1c: Force Per Unit Span vs. Angle of Attack");
legend('Location','best')

subplot(2,2,4)
plot(alpha,ld);
hold on
grid on
xlabel("Angle of Attack (deg)");
ylabel("Lift/Drag Ratio");
title("Fig 1d: Lift/Drag Ratio vs. Angle of Attack");

hold off
fprintf("1 figure has been outputted.");

fprintf('\nKyra Bryan''s AEM413 Project: 4.17 script complete. --------------------------------------------------\n\n');