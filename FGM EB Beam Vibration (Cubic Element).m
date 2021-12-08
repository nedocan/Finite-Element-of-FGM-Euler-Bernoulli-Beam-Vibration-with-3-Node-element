clc
clear

syms x L c1 c2 c3 c4 c5 c6 c7 c8 c9...
    W1 W2 W3 Teta1 Teta2 Teta3 U1 U2 U3...
    A11 B11 D11 I0 I1 I2 z real

w=c1+c2*x+c3*x^2+c4*x^3+c5*x^4+c6*x^5;
u=c7+c8*x+c9*x^2;

u1=subs(u,0);
w1=subs(w,0);
teta1=subs(diff(w,x),0);
u2=subs(u,L/2);
w2=subs(w,L/2);
teta2=subs(diff(w,x),L/2);
u3=subs(u,L);
w3=subs(w,L);
teta3=subs(diff(w,x),L);

eqns=[u1
      w1
      teta1
      u2
      w2
      teta2
      u3
      w3
      teta3];
  
q=[U1 W1 Teta1 U2 W2 Teta2 U3 W3 Teta3];
vars=[c1 c2 c3 c4 c5 c6 c7 c8 c9];
[coef] = equationsToMatrix(eqns,vars);
c=inv(coef)*q';

w=collect(c(1)+c(2)*x+c(3)*x^2+c(4)*x^3+c(5)*x^4+c(6)*x^5,q);
u=collect(c(7)+c(8)*x+c(9)*x^2,q);

[N_u_local,terms]=coeffs(u,q);
[N_w_local,terms]=coeffs(w,q);
N_w=[0 N_w_local(1) N_w_local(2) 0 N_w_local(3) N_w_local(4) 0 N_w_local(5) N_w_local(6)];
N_u=[N_u_local(1) 0 0 N_u_local(2) 0 0 N_u_local(3) 0 0];

K_u=A11*int(diff(N_u,x)'*diff(N_u,x),x,[0,L]);
K_w=D11*int(diff(N_w,x,x)'*diff(N_w,x,x),x,[0,L]);
K_u_w=-B11*int(diff(N_u,x)'*diff(N_w,x,x)+diff(N_w,x,x)'*diff(N_u,x),x,[0,L]);
K=K_u+K_w+K_u_w;

M_u=I0*int(N_u'*N_u,x,[0,L]);
M_w=I0*int(N_w'*N_w,x,[0,L])+I2*int(diff(N_w,x)'*diff(N_w,x),x,[0,L]);
M_u_w=-I1*int(N_u'*diff(N_w,x)+diff(N_w,x)'*N_u,x,[0,L]);
M=M_u+M_w+M_u_w;

%%%% Geometry %%%%%
h=  0.1; %height of beam
b=   0.05; %width of beam
l=   2;%length of beam  

%%% Shear Considerations %%%
f=0 ;   %Euler Bernoulli Beam

%%%% Material Properties %%%%%
k=2;  %volume fraction

E_metallic=210*10^9;
rho_metallic=7800;
nu_metallic=0.3;

E_ceramic=380*10^9;
rho_ceramic=3800;
nu_ceramic=0.14;

%%%%%  Case; upper face is ceramic; lower face is metallic %%%%%
E_effective=(E_ceramic-E_metallic)*(z/h+0.5)^k+E_metallic;
rho_effective=(rho_ceramic-rho_metallic)*(z/h+0.5)^k+rho_metallic;
nu_effective=(nu_ceramic-nu_metallic)*(z/h+0.5)^k+nu_metallic;

%%% Inertial Properties %%%
I0_numeric=double(int(rho_effective,[-h/2 h/2])*b);
I1_numeric=double(int(z*rho_effective,[-h/2 h/2])*b);
I2_numeric=double(int(z^2*rho_effective,[-h/2 h/2])*b);
I=[I0_numeric,I1_numeric,I2_numeric] ;%% Ineartias

%%%Coefficients of Stress Resultants%%%
A11_numeric=double(int(E_effective,[-h/2 h/2])*b);
B11_numeric=double(int(z*E_effective,[-h/2 h/2])*b);
D11_numeric=double(int(z^2*E_effective,[-h/2 h/2])*b);
R=[A11_numeric,B11_numeric,D11_numeric]; %% Coefficinet of Stress Resultants

nel=5;  %%Number of Divided Elements

mm=double(subs(M,{A11,B11,D11,I0,I1,I2,L},{R(1),R(2),R(3),I(1),I(2),I(3),l/nel}));
kk=double(subs(K,{A11,B11,D11,I0,I1,I2,L},{R(1),R(2),R(3),I(1),I(2),I(3),l/nel}));

n_nodes=2*nel+1;  %%Node Number
Dof=3;  %%Degree of Freedom of a Node

%%% Boundary Conditions %%%
BCs=[1,2,3,n_nodes*Dof-2,n_nodes*Dof-1,n_nodes*Dof]; % C-C
%BCs=[1,2,n_nodes*Dof-2,n_nodes*Dof-1]; % S-S
%BCs=[1,2,3,n_nodes*Dof-2,n_nodes*Dof-1]; %C-S
%BCs=[1,2,3]; % C-F



           
Kglob=zeros(Dof*n_nodes,Dof*n_nodes)  ;
Mglob=zeros(Dof*n_nodes,Dof*n_nodes)  ;

for i=1:nel
    top(i,:)=[6*i-5 6*i-4 6*i-3 6*i-2 6*i-1 6*i 6*i+1 6*i+2 6*i+3]; %% Topology Matrix
end

for n=1:nel
    for i=1:9
        for j=1:9
            Kglob(top(n,i),top(n,j))=...
                Kglob(top(n,i),top(n,j))+kk(i,j);
        end
    end
end

for n=1:nel
    for i=1:9
        for j=1:9
            Mglob(top(n,i),top(n,j))=...
                Mglob(top(n,i),top(n,j))+mm(i,j);
        end
    end
end

activeDof=setdiff(1:Dof*n_nodes',BCs);

Kglobactive=Kglob(activeDof,activeDof);
Mglobactive=Mglob(activeDof,activeDof);

%System_Matrix=(inv(Mglobactive)*Kglobactive);
[vectorc,freqc]=eig(Kglobactive,Mglobactive);
nat_freq=freqc;
freqc=diag(freqc);
freqc=sort(sqrt(freqc))

