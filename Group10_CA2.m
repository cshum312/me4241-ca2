% CA2 Submission-Group10
options = optimset('Display','off');

global S mass W R P_SL T_SL rho_SL h_tropo c
S = 17.1;
mass = 1248.5;
W = mass * 9.81;
c = 1.74;


R = 287;             
P_SL = 101325;         
h_tropo = 11000;           
T_SL = 15 + 273.15;  
rho_SL = P_SL/(R*T_SL);

h_array = [];
V_array = []; 
de_array = [];

alpha_max = 25*pi/180;
alpha = alpha_max;
d_alpha = alpha_max/100;
x0 = [0,100,0];

while alpha >= d_alpha
    
    fun_1 = @(x)solve_flight_env(x,alpha);
    x = fsolve(fun_1,x0, options);
    alpha = alpha - d_alpha;
    x0 = x;
    
    if (rho(x(1)) > 1.2) || (rho(x(1)) < 0.1)
%         disp("rho");
%         disp(rho(x(1)));
        continue
        
    end

    if (rad2deg(x(3)) > 30) || (rad2deg(x(3)) < -30) 
%         disp("elevator");
%         disp(rad2deg(x(3)));
        continue
    end

    h_array(end+1) = x(1)/1000;
    V_array(end+1) = x(2);
    de_array(end+1) = x(3);

end

disp("2. Refer to figure 1")
figure(1)
plot(V_array, h_array);
title('Level Flight Envelope (h-V Plot)');
xlabel('Velocity (m/s)');
ylabel('Altitude (km)');
grid on

for n=4:length(h_array)   
    c = 1.74;
    S = 17.1;
    mass = 1248.5;
    W = mass * 9.81;
    
    Iy = 4067.5;
    g = 9.81;
    h = h_array(n);                  
    
    Rho = rho(h);                                                           
    de = de_array(n);                 
    
    syms alphaT VT qT thetaT
    
    %Non-linear EOM
    eqn1 = qT -(T(VT,h,alphaT)*sin(alphaT)+0.5*Rho*(VT/cos(alphaT))^2*S*CL(alphaT,qT,c,VT,de))/(mass*VT)+(g/VT)*cos(thetaT-alphaT);
    eqn2 = (T(VT,h,alphaT)*cos(alphaT)-0.5*Rho*(VT/cos(alphaT))^2*S*CD(alphaT,qT,c,VT,de))/(mass)-g*sin(thetaT-alphaT);
    eqn3 = 0.5*Rho*(VT/cos(alphaT))^2*S*c*Cm(alphaT,qT,c,VT,de)/Iy;
    eqn4 = qT; 
    
    
    J(alphaT,VT,qT,thetaT) = jacobian([eqn1,eqn2,eqn3,eqn4],[alphaT,VT,qT,thetaT]);
    sol = vpasolve([eqn1==0, eqn2==0, eqn3==0, eqn4==0],[alphaT VT qT thetaT],[0.113;123;0;0.113]);
    
    %aoa,u,q,theta at trim/fixed point
    alphaT = sol.alphaT;                  
    VT = sol.VT;                       
    qT = sol.qT;                       
    thetaT = sol.thetaT;    
    A = J(alphaT,VT,qT,thetaT);

    %short period and phugoid eigenvalues
    eigenvalues = eig(A);            
    eigenvalues_array(n,:) = eigenvalues;    

end

eigenvalues_array;
disp("3. Refer to Figure 2")
figure(2)
plot(eigenvalues_array(:,1),"r.")
hold on
plot(eigenvalues_array(:,3),"b.")
hold on
plot(eigenvalues_array(:,2),"r.")
hold on
plot(eigenvalues_array(:,4),"b.")
hold off
grid on
title('Eigenvalues Plot');
xlabel('real');
ylabel('imaginary');
legend("Phugoid Mode","Short Period Mode")

disp("Comment on the stability of the drone:")
disp("As shown in the eigenvalues plot, the phugoid mode has multiple points (four points) with real positive eig values, this is also where the UAV is dynamically unstable; does not recover from pertubations.")
disp("Due to these points, a flight controller can be considered as it allows the UAV to reduce the damping effects and improves stability.")

function F = solve_flight_env(x,alpha)  

    global S W c
    %q=0
    %x(1) is h, x(2) is V_T, x(3) is de
    CL = 6.44*alpha + 0.355*x(3);
    CD = 0.03 + 0.05*CL^2;
    Cm = 0.05 - 0.683*alpha - 0.923*x(3);
    Power = 59164 + 64304*sqrt(rho(x(1))) - 15312*rho(x(1));
    
    T = Power/x(2);                          %Thrust = Power/velocity
    L =  0.5*rho(x(1))*x(2)^2*S*CL;
    D = 0.5*rho(x(1))*x(2)^2*S*CD;
    m = 0.5*rho(x(1))*x(2)^2*S*c*Cm;
    
    F(1) = T*cos(alpha)-D;
    F(2) = L+T*sin(alpha)-W;
    F(3) = m;

end

function atmo = atmo(h)
    t = 11000;        
    T_sl = 15;        
   
    c = 6.51/1000; 
    R = 287.058; 
    g = 9.81; 
    x = g/(R*c); 
    
    %SL conditions
    T_sl = T_sl + 273.15; 
    P_sl = 101325; 
    R_sl = P_sl/(R*T_sl);
    
    %Tropopause
    T_st = T_sl - c * t; 
    P_st = P_sl * (T_st/T_sl)^x; 
    R_st = R_sl * (T_st/T_sl)^(x-1); 
    
    if h <= t 
        T = T_sl - c * h; 
        P = P_sl * (T/T_sl)^x; 
        Rho = R_sl * (T/T_sl)^(x-1); 
        a = sqrt(1.4*R*T);
        
    else    
        T = T_st;
        P = P_st*exp(-g/(R*T_st)*(h-t));
        Rho = R_st*exp(-g/(R*T_st)*(h-t));
        a = sqrt(1.4*R*T);
    end
       atmo = [T P Rho a];
end

function Rho = rho(h)
    atmo_h = atmo(h);
    Rho = atmo_h(3);
end

function CL = CL(aoa,q,c,u,de)
    CL = 6.44*aoa+3.8*(q*c)/(2*u)+0.355*de;
end

function CD = CD(aoa,q,c,u,de)
    CD =0.03 + 0.05*(CL(aoa,q,c,u,de))^2;
end

function Cm = Cm(aoa,q,c,u,de)
    Cm =0.05 - 0.683*aoa - 9.96*(q*c)/(2*u) - 0.923*de;
end

function T = T(V,h,aoa)
    Rho = rho(h);
    T = (59164 + 64304*sqrt(Rho) - 15312*(Rho))/(V*cos(aoa));
end