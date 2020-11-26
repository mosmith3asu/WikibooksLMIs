[Zf,Pf,gammaf]=function(A,B,C,options)

% Define sizes of matricies
p = size(C,1);
n= size(A,1);
m = size(B,2);

% Define variables
Z = sdpvar(n,p);
P = sdpvar(n);
gamma = sdpvar(1);

% Define matricies
m11=[-P(A+B*F)+(A+B*F)'*P+Z*C+C'*Z'];
m12 = [-Z];
m13 = [X*B];
m22 = [-gamma*eye(n)];
m23 = [zeros(n)];
m33 = [-gamma*eye(n)];

Mat1 = [P*A+A'+Z*C+C'*Z'];
Mat2 = [
    m11 m12 m3;
    m12' m21 m23;
    m13' m23' m33];

% Compile constraints
Constraints = [
    gamma>eps*eye(1);
    Mat1 < eps*eye(size(Mat1));
    Mat2< eps*eye(size(Mat2))];

objective = gamma;
optimize(Constraints, Obective, options);

gammaf = value(gamma);
Pf = value(P);
Zf = value(Z);

end