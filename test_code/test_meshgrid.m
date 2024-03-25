% 用来测试三维画图

num = 150;
a = linspace(0, 2*pi, num);
b = linspace(-0.5*pi, 0.5*pi, num);
[a, b] = meshgrid(a, b);
r = 5;
c = sqrt(abs(a - pi))*1.5;
X = r*cos(b).*sin(a).*c;
Z = -r*cos(b).*cos(a).*c;
Y = r.*sin(b)*0.75;
for i = 1:num
    for j =1:num
        C(i,j,1) = Z(i, j)/(max(max(Z))-min(min(Z))) + abs(min(min(Z)))/(max(max(Z))-min(min(Z)));
        C(i,j,2) = 0;
        C(i,j,3) = 0;
    end
end
s = surf(X, Y, Z, C);
axis equal
axis tight
s.FaceAlpha = 0.9;
s.EdgeColor = 'none';
s.FaceColor = 'interp';


%%
figure(1)
x=-2:0.1:2;   y=-3:0.1:3; 
[X,Y]=meshgrid(x,y);
Z=X.*exp(-X.^2-Y.^2);
mesh(X,Y,Z);