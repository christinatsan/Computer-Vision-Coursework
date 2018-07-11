

%% Generate 3D calibration pattern: 
%% Pw holds 32 points on two surfaces (Xw = 1 and Yw = 1) of a cube 
%%Values are measured in meters.
%% There are 4x4 uniformly distributed points on each surface.

cnt = 1;

%% plane : Xw = 1

for i=0.2:0.2:0.8,
 for j=0.2:0.2:0.8,
   Pw(cnt,:) = [1 i j];
   cnt = cnt + 1;
 end
end

%% plane : Yw = 1

for i=0.2:0.2:0.8,
 for j=0.2:0.2:0.8,
   Pw(cnt,:) = [i 1 j];
   cnt = cnt + 1;
 end
end

N = cnt;

%plot3(Pw(:,1), Pw(:,2), Pw(:,3), '+');

%% Virtual camera model 

 %% Extrinsic parameters : R = RaRbRr

gamma = 40.0*pi/180.0;
Rr = [ [cos(gamma) -sin(gamma) 0];
       [sin(gamma) cos(gamma)  0];
       [  0          0         1]; ];

beta = 0.0*pi/180.0;

Rb = [ [cos(beta) 0 -sin(beta)];
       [0         1       0];
       [sin(beta) 0  cos(beta)]; ];

alpha = -120.0*pi/180.0;
Ra = [ [1      0                0];
       [0   cos(alpha)  -sin(alpha)];
       [0   sin(alpha)   cos(alpha)]; ];

R = Ra*Rb*Rr;

T = [0 0 4]';

%% Intrinsic parameters

f = 0.016;
Ox = 256;
Oy = 256;

Sx = 0.0088/512.0;
Sy = 0.0066/512.0;

Fx = f/Sx;
Fy = f/Sy;

%% asr is the aspect ratio
asr = Fx/Fy;

%% Generate Image coordinates

%% surface Xw = 1
cnt = 1;
for cnt = 1:1:16
   Pc(cnt,:) = (R*Pw(cnt,:)' + T)';
   p(cnt,:)  = [(Ox - Fx*Pc(cnt,1)/Pc(cnt,3)) (Oy - Fy*Pc(cnt,2)/Pc(cnt,3))];
end
plot(p(:,1), p(:,2), 'r+');
axis([0 512 0 512]);
hold;

%% surface Yw = 1
for cnt = 17:1:32,
   Pc(cnt,:) = (R*Pw(cnt,:)' + T)';
   p(cnt,:)  = [(Ox - Fx*Pc(cnt,1)/Pc(cnt,3)) (Oy - Fy*Pc(cnt,2)/Pc(cnt,3))];
end
plot(p(17:32,1), p(17:32,2), 'g+');
grid;
figure;

plot3(Pc(:,1), Pc(:,2), Pc(:,3), '+');
grid;
figure;

%%

X = Pw(:,1);
Y = Pw(:,2);
Z = Pw(:,3);

u_values = p(:,1)-256;
v_values = p(:,2)-256;

% Build matrix A
A = [u_values.*X, u_values.*Y, u_values.*Z, u_values, -v_values.*X, -v_values.*Y, -v_values.*Z, -v_values];

% Get solution v
[U, D, V] = svd(A,0);
%V = V.';
v = V(:,end);

% Determine scale
abs_gamma = sqrt(v(1)^2 + v(2)^2 + v(3)^2);

% Determine aspect ratio
a = sqrt(v(5)^2+v(6)^2+v(7)^2)/sqrt(v(1)^2+v(2)^2+v(3)^2);

% Determine 1st two rows of R

R = zeros(3,3);
R(1,1) = v(5)/(abs_gamma*a);
R(1,2) = v(6)/(abs_gamma*a);
R(1,3) = v(7)/(abs_gamma*a);
R(2,1) = v(1)/abs_gamma;
R(2,2) = v(2)/abs_gamma;
R(2,3) = v(3)/abs_gamma;

% Determine 1st two components of T

T = zeros(3,1);
T(1) = v(8)/(abs_gamma*a);
T(2) = v(4)/abs_gamma;

% Determine sign s

Xc = R(1,1)*Pw(1,1) + R(1,2)*Pw(1,2) + R(1,3)*Pw(1,3) + T(1);
xc = p(1,1);

if Xc*xc < 0
    s = -1;
else
    s = 1;
end

% Recalculate using sign s

R(1,1) = s*R(1,1);
R(1,2) = s*R(1,2);
R(1,3) = s*R(1,3);
R(2,1) = s*R(2,1);
R(2,2) = s*R(2,2);
R(2,3) = s*R(2,3);
T(1) = s*T(1);
T(2) = s*T(2);

% Compute 3rd row of R

R3_transpose = cross(R(1,:).',R(2,:).');
R(3,1) = R3_transpose(1);
R(3,2) = R3_transpose(2);
R(3,3) = R3_transpose(3);

% Enforce orthogonality constraint on R

[U,D,V2] = svd(R);
R = U*V2.';

% Solve for Tz and fx and fy

A2 = [u_values, R(1,1).*X + R(1,2).*Y + R(1,3).*Z + T(1)];
B2 = -1*u_values.*(R(3,1).*X + R(3,2).*Y + R(3,3).*Z);

sol = (pinv(A2.'*A2))*A2.'*B2;

T(3) = sol(1);
fx = sol(2);
fy = fx/a;

ox = 256;
oy = 256;

% Find Projective Matrix

Mint = [-fx 0 ox; 0 -fy oy; 0 0 1];
Mext = [R(1,1) R(1,2) R(1,3) T(1); R(2,1) R(2,2) R(2,3) T(2); R(3,1) R(3,2) R(3,3) T(3)];
M = Mint*Mext;

% Reproject
worldCoords = [X Y Z];
points4D = [worldCoords ones(size(worldCoords,1),1) ];

reproj = M*points4D.';
newMatrix = zeros(2,size(reproj,2));
newMatrix(1,:) = reproj(1,:)./reproj(3,:);
newMatrix(2,:) = reproj(2,:)./reproj(3,:);
newMatrix = newMatrix.';

plot(newMatrix(:,1), newMatrix(:,2), 'b+');

grid;
axis([0 512 0 512]);

