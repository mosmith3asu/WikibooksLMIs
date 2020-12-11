clc
clear all
% ==============================
% Define System
Ad=[ 0 1 0;-30 0 20;2000 0 -2000];
Bd=[0;0.1;0];
Cd=[0 0 1];
Dd = 0;


n=size(Ad,1); % number of states
m=size(Bd,2);
p=size(Cd,1);
% ==============================
% Begin Optimization

gamma=sdpvar(1);
P=sdpvar(n);
Z=sdpvar(n,p);
F=sdpvar(m,n);
X=sdpvar(n,n);
N11=sdpvar(n,n);

MAT1=[
    Ad'*P*Ad-P-Ad'*Z*Cd-Cd'*Z'*Ad Cd'*Z';
    (Cd'*Z')' -P];

M11=N11; M12=(Ad+Bd*F)'*Z; M13 =X*Bd; M14 = Cd'*Z';
M21 = M12'; M22 = -gamma*ones(p,1); M23 = zeros(p); M24 = Z';
M31 = M13'; M32 = M23'; M33 = -gamma*ones(m); M34 = zeros(m,n);
M41 = M14'; M42 = M24'; M43 =zeros(n,1); M44 = -P;

MAT2=[
    M11 M12 M13 M14;
    M21 M22 M23 M24;
    M31 M32 M33 M34;
    M41 M42 M43 M44];
% MAT2=[
%     N11 (Ad+Bd*F)'*Z X*Bd Cd'*Z';
%     ((Ad+Bd*F)'*Z)' -gamma*ones(p,1) zeros(p) Z';
%     (X*Bd)' zeros(m) -gamma*ones(m) zeros(m);
%     (Cd'*Z')' Z zeros(n) -P];

F=[F; P>=0];
F=[F; MAT1<=eps*eye(size(MAT1))];
F=[F; MAT2<=eps*eye(size(MAT2))];

sol = optimize(F);

% ==============================
% Report
feasable= sol.problem;
if feasable==0
    disp('Strongly Stabalizable')
else
    disp('NOT Strongly Stabalizable')
end