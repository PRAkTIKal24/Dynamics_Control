%% Linearizing the cart pole system
clc
clear
%Defining required parameters as syms
syms x1 x2 x3 x4 u m M g l

a = 1/(m+M);
%State vector definition
x = [x1;x2;x3;x4];

%Input Matrices from non-linear state space equation
Atemp = [x2;
     ((-m*a*g*sin(2*x3)/2)+(a*l*x4^2*sin(x3)*2*m))/(2-m*a*cos(x3)^2);
     x4;
     (g*sin(x3)-m*l*a*x4^2*sin(2*x3)/2)/((2*l)-(m*l*a*cos(x3)^2))];
 
btemp = [0;
     (2*a*u)/(2-m*a*cos(x3)^2);
     0;
     -a*cos(x3)*u/(2*l-m*l*a*cos(x3)^2);];

%Calculating matrices corresponding to the linear state space form
A1 = jacobian(Atemp,x);
B1 = jacobian(btemp,u);

A1 = subs(A1, {x1,x2,x3,x4}, {0,0,0,0});
B1 = subs(B1, {x1,x2,x3,x4,u}, {0,0,0,0,0});

A1 = simplify(expand(A1))
B1 = simplify(expand(B1))

%So the linearized model is x_dot = A1x + B1u


%% Designing state controller

%Defining gain and pole symbols
syms k1 k2 k3 k4 S;

%Gain, pole vector definition
K = [k1 k2 k3 k4];

j = sqrt(-1); %complex number def
poles = [-1 -2 -1+j -1-j]; %Given poles

%Identity matrix
I = eye(4);

%Calculating Det(A-bK-sI), where s are the poles
Determ = det(A1-B1*K-S*I);

Determ = simplify(expand(Determ))

%Substitute the pole values to form 4 equations, solving which will give
%the gain values for the given poles
for i=1:3
    Eqn(i) = subs(Determ, {S,M,m,l,g}, {poles(i),1,0.1,1,10});
    Eqn(i) == 0; %Equation definition
end

%Solve the equations to get gain values
sol = solve(Eqn,K);
K1 = simplify(expand(sol.k1));
K2 = simplify(expand(sol.k2));
K3 = simplify(expand(sol.k3));
K4 = simplify(expand(sol.k4));
%Print out gain matrix
Gain = [K1 K2 K3 K4]

%% Cart pole Animation based on above state controller
%Inputting initial position, orientaton values from user
prompt = 'Enter initial position value followed by the initial anglular orientation of point mass wrt verticle: ';
X = input(prompt)
Theta = input(prompt)

%The controlled model is z_dot = (A1-B1*K)*z
%Initializing new variables to not interfere with previous sections
syms t z1(t) z2(t) z3(t) z4(t)

%Defining z,z_dot matrices
z = [z1(t);z2(t);z3(t);z4(t)];
z_dot = diff(z,t);
z2(t) = diff(z1,t);
z4(t) = diff(z4,t);

%Defining state controlled model eqn,
Mm = A1-B1*Gain; %Just to subs in the values of masses, length and 'g'
Mm = simplify(expand(Mm));
Mm = subs(Mm, {m,M,l,g}, {0.1,1,1,10});

eqn = z_dot - Mm*z;
eqn == 0;

%Initial conditions
cond1 = z1(0) == X;
cond2 = z3(0) == Theta; 
conds = [cond1; cond2];

%Solving the differential equations for z1,z2,z3,z4
Zz = dsolve(eqn, conds);

%Plot the generated Z values to simulate the movement of the pole under
%state feedback control

%Solve the equation for 'time' time steps
time = 10;

f = figure;
for i=0:time
    plot(subs(Zz.z1, {t},{i}),i, '-o');
    hold on;
    plot([subs(Zz.z1, {t},{i})+sin(subs(Zz.z3, {t},{i}))],[0.5+cos(subs(Zz.z3, {t},{i}))], '-X');
    rectangle('Position', [subs(Zz.z1, {t},{i})-1, subs(Zz.z1, {t},{i})+1, 0.5, 0.5], 'EdgeColor', 'b');
    hold off;
end
    
    
