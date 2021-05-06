clc
clear 
close all

%% Initial Comments 

%{
Kyra Bryan
AEM 313-001 
Project 2: Airfoil Study

TAT
Calculate angle of zero lift for 4 airfoils
Calculate moment about c/4 for 4 airfoils
Calculate moment about LE for 4 airfoils
Calculate center of pressure for 4 airfoils

XFOIL
Inviscid solutions
Viscous solutions (Re=10^6)
Variation of Re above and below prescribed value (one airfoil)
Variation of thickness (two airfoils)
%} 

%% TAT

fprintf('*****TAT*****\n\n');

airfoil = ["2321","2421","4321","4412"];
p = [0.3,0.4,0.3,0.4];
m = [0.02,0.02,0.04,0.04];
c=1; %to keep the (x/c) consistent in the below equations

for i = 1:4 
    
    %given eqn for mean camber line
    syms x theta
    z1 = ((m(i)*c)/p(i)^2)*(2*p(i)*(x/c)-(x/c)^2);
    dz1_cart = diff(z1,x); %here, x=x/c
    z2 = ((m(i)*c)/(1-p(i))^2)*((1-2*p(i))+2*p(i)*(x/c)-(x/c)^2);
    dz2_cart = diff(z2,x); %here, x=x/c
    
    %transform cartesian to polar with x = (c/2)*(1-cos(theta)) (from section 4.8)
    theta_split = acos(1-(2*p(i)/c)); %dz1 for 0<=theta<=theta_split, dz2 for theta_split<=theta<=pi
    dz1 = subs(dz1_cart,x,(c/2)*(1-cos(theta)));
    dz2 = subs(dz2_cart,x,(c/2)*(1-cos(theta)));
    
    %Fourier coefficients, alpha0, cmc4
    syms alpha
    F_A1_1 = dz1*cos(theta);
    F_A1_2 = dz2*cos(theta);
    A1 = (2/pi)*int(F_A1_1,theta,0,theta_split) + (2/pi)*int(F_A1_2,theta,theta_split,pi); %eqn 4.51
    F_A2_1 = dz1*cos(2*theta);
    F_A2_2 = dz2*cos(2*theta);
    A2 = (2/pi)*int(F_A2_1,theta,0,theta_split) + (2/pi)*int(F_A2_2,theta,theta_split,pi); %eqn 4.51
    F_alpha0_1 = dz1*(cos(theta)-1);
    F_alpha0_2 = dz2*(cos(theta)-1);
    alpha0 = -(1/pi)*int(F_alpha0_1,theta,0,theta_split) + -(1/pi)*int(F_alpha0_2,theta,theta_split,pi); %eqn 4.61
    cmc4 = (pi/4)*(A2-A1); %eqn 4.64
    
    %output values into command window
    fprintf('For the NACA %s airfoil:\n',airfoil(i));
    fprintf('Angle of zero lift = %.6f rad = %.6f deg\n',alpha0,rad2deg(double(alpha0)));
    fprintf('Moment coefficient about c/4 = %.6f\n',cmc4);
    
    %AOA-dependent values to plot variation with AOA
    alpha = -6:2:30;
    cl = 2*pi*(deg2rad(alpha)-double(alpha0)); %eqn 4.57    
    cmle = -((cl/4)+(pi/4)*(A1-A2)); %eqn 4.63
    xcp = (c/4)*(1+(pi./cl)*(A1-A2)); %eqn 4.66
    
    %plot cmle
    fig = figure();
    plot(alpha,cmle);
    grid on    
    title("TAT - c_m_,_l_e vs. AOA - NACA " + airfoil(i));
    xlabel('Angle of Attack (deg)');
    xticks(alpha);
    xlim([-8 32]);
    ylabel('Coefficient of Moment about LE');
    saveas(fig,airfoil(i)+'_cmle.png')
    
    %plot xcp
    fig = figure();
    plot(alpha,xcp);
    grid on
    title("TAT - x_c_p vs. AOA - NACA " + airfoil(i));
    xlabel('Angle of Attack (deg)');
    xticks(alpha);
    xlim([-8 32]);
    ylabel('Center of Pressure');
    saveas(fig,airfoil(i)+'_xcp.png')
    
    fprintf('\n');

end 

fprintf('TAT plots have been generated and saved to the working directory.\n\n');

%% XFOIL

fprintf('*****XFOIL*****\n\n');

%data from XFOIL
inviscid = [-0.522400000000000,-0.263900000000000,-0.00520000000000000,0.253600000000000,0.512100000000000,0.770000000000000,1.02690000000000,1.28260000000000,1.53670000000000,1.78890000000000,2.03900000000000,2.28660000000000,2.53140000000000,2.77310000000000,3.01150000000000,3.24620000000000,3.47700000000000,3.70350000000000,3.92550000000000;
            -0.503600000000000,-0.245100000000000,0.0137000000000000,0.272500000000000,0.530900000000000,0.788700000000000,1.04560000000000,1.30110000000000,1.55510000000000,1.80720000000000,2.05710000000000,2.30450000000000,2.54910000000000,2.79060000000000,3.02870000000000,3.26320000000000,3.49360000000000,3.71980000000000,3.94150000000000;
            -0.270700000000000,-0.0115000000000000,0.247700000000000,0.506600000000000,0.764900000000000,1.02220000000000,1.27830000000000,1.53290000000000,1.78550000000000,2.03610000000000,2.28410000000000,2.52940000000000,2.77150000000000,3.01040000000000,3.24550000000000,3.47670000000000,3.70370000000000,3.92620000000000,4.14390000000000;
            -0.216300000000000,0.0259000000000000,0.268100000000000,0.510000000000000,0.751200000000000,0.991500000000000,1.23070000000000,1.46830000000000,1.70410000000000,1.93790000000000,2.16930000000000,2.39800000000000,2.62390000000000,2.84650000000000,3.06580000000000,3.28120000000000,3.49270000000000,3.70000000000000,3.90270000000000];
viscous = [-0.460700000000000,-0.233200000000000,-0.00340000000000000,0.228300000000000,0.457100000000000,0.686900000000000,0.908500000000000,1.11950000000000,1.30220000000000,1.46320000000000,1.61210000000000,1.72800000000000,1.80880000000000,1.83500000000000,1.79650000000000,1.71250000000000,1.64460000000000,NaN,1.58480000000000;
           -0.447600000000000,-0.218600000000000,0.0130000000000000,0.245800000000000,0.475800000000000,0.699600000000000,0.916800000000000,1.12090000000000,1.29560000000000,1.45740000000000,1.60020000000000,1.71780000000000,1.80110000000000,1.83500000000000,1.82060000000000,1.76950000000000,1.69920000000000,1.65690000000000,1.61510000000000;
           -0.244700000000000,-0.0157000000000000,0.216000000000000,0.449200000000000,0.680500000000000,0.906900000000000,1.12780000000000,1.32950000000000,1.51420000000000,1.67580000000000,1.80290000000000,1.88600000000000,1.91880000000000,1.89190000000000,1.79590000000000,1.71460000000000,1.67760000000000,1.66200000000000,NaN;
           -0.207900000000000,0.0222000000000000,0.253000000000000,0.484400000000000,0.714500000000000,0.934700000000000,NaN,NaN,NaN,1.70810000000000,1.83760000000000,1.93680000000000,1.99200000000000,2.00970000000000,1.95870000000000,1.86080000000000,1.77440000000000,1.68440000000000,NaN];

rey1 = [10,100,1000,1e4,1e5,1e6,10e6,1e8,1e9,1e10];
cl_re1 = [0.2889,0.1935,0.0183,-0.0421,1.0080,0.8223,0.9085,0.9318,0.9493,0.9674];
rey2 = [5e5,1e6,5e6,1e7,5e7,1e8,5e8];
cl_re2 = [0.8181,0.8223,0.8951,0.9085,0.9265,0.9318,0.9440];
rey3 = [7e6,8e6,9e6,10e6,11e6,12e6,13e6];
cl_re3 = [0.9021,0.9049,0.9074,0.9085,0.9107,0.9120,0.9125]; 

%plot Re plots
fprintf('Generating Re variation plots...\n\n');

fig = figure();
semilogx(rey1,cl_re1);
grid on
title("Re Variation - Log Scale - NACA " + airfoil(1));
xlabel('Reynolds number');
ylabel('Coefficient of lift');
saveas(fig,airfoil(1)+'_Re1.png')

fig = figure();
semilogx(rey2,cl_re2);
grid on
title("Re Variation - Log Scale, Re>1e5 only - NACA " + airfoil(1));
xlabel('Reynolds number');
ylabel('Coefficient of lift');
saveas(fig,airfoil(1)+'_Re2.png')

fig = figure();
plot(rey3,cl_re3);
grid on
title("Re Variation - Even Scale - NACA " + airfoil(1));
xlabel('Reynolds number');
ylabel('Coefficient of lift');
saveas(fig,airfoil(1)+'_Re3.png')

fprintf('Cp and thickness plots were generated in XFOIL.\n\n');
fprintf('XFOIL plots have been generated and saved to the working directory.\n\n');

%% Comparisons

fprintf('*****Comparisons*****\n\n');

fprintf('Generating comparison plots...\n\n');

for i = 1:4
    
    %plot TAT vs XFOIL
    fig = figure();
    plot(alpha,cl,'DisplayName','TAT Calculation','Color', 'r');
    hold on
    grid on
    plot(alpha,inviscid(i,:),'DisplayName','XFOIL Inviscid Solution','Color', 'b');
    legend('Location','best')
    title("TAT vs. XFOIL - c_l vs. AOA - NACA " + airfoil(i));
    xlabel('Angle of Attack (deg)');
    xticks(alpha);
    xlim([-8 32]);
    ylabel('Coefficient of Lift');
    saveas(fig,airfoil(i)+'_TATXFOILCompare.png')
    hold off
    
    %plot inv and vis
    fig = figure();
    plot(alpha,inviscid(i,:),'DisplayName','Inviscid Solution','Color', 'b');
    hold on
    grid on
    plot(alpha,viscous(i,:),'DisplayName','Viscous Solution','Color', 'r');
    legend('Location','best')
    title("XFOIL - Inviscid vs. Viscous - c_l vs. AOA - NACA " + airfoil(i));
    xlabel('Angle of Attack (deg)');
    xticks(alpha);
    xlim([-8 32]);
    ylabel('Coefficient of Lift');
    saveas(fig,airfoil(i)+'_InvVisCompare.png')
    hold off
    
end 

fprintf('Comparison plots have been generated and saved to the working directory.\n\n');

%% Final Output

%close all %uncomment to close all pop-up windows of figures

fprintf('Script completed.\n\n');
