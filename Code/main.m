clc
clear all
close all
img = imread('000852.jpg');
model_len=[160,95,95,65,65,60];
part = 1; seq = 7;
[m,n]=size(img);
x = [0:10:n]; y = [0:10:m]; theta = pi/2; s = 1;
lx = length(x); ly = length(y);
min_e = inf;
lF = ReadStickmenAnnotationTxt('buffy_s5e2_sticks.txt');
tic
for i = 1:lx
    for j = 1:ly
        z = [x(i),y(j),theta,s];
        cost = match_energy_cost(z,part,seq,lF);
        if cost < min_e
            min_e = cost;
            best_i = i; best_j = j;
        end
    end
end

x = x(best_i); y = y(best_j);
theta = [-pi/2:pi/20:pi/2]; s = [0.2:0.2:0.8 1:0.2:3];
lt = length(theta); ls = length(s);
min_e = inf;
for p = 1:lt
    for q = 1:ls
        z = [x,y,theta(p),s(q)];
        cost = match_energy_cost(z,part,seq,lF);
        if cost < min_e
            min_e = cost;
            best_p = p; best_q = q;
        end
    end
end
z = [x,y,theta(best_p),s(best_q)];

l1 = [z(1)-0.5*z(4)*model_len(part)*sin(z(3));
    z(2)+0.5*z(4)*model_len(part)*cos(z(3));
    z(1)+0.5*z(4)*model_len(part)*sin(z(3));
    z(2)-0.5*z(4)*model_len(part)*cos(z(3))];
[l2,joint2] = findl(z,1,2,seq,m,n,[]);
[l3,joint3] = findl(z,1,3,seq,m,n,[]);
l6 = findl(z,1,6,seq,m,n,[]);
l4 = findl(z,2,4,seq,m,n,joint2);
l5 = findl(z,3,5,seq,m,n,joint3);
L = [l1,l2,l3,l4,l5,l6];
DrawStickman(L,img)
toc