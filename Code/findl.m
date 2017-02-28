function [l,joint] = findl(z,z_part,w_part,seq,m,n,joint1)
% z is the parent(root) location, for part, 2 = left upper arm, 
% 3 = right upper arm, 6 = head
lF = ReadStickmenAnnotationTxt('buffy_s5e2_sticks.txt');
L = 130;
model_len=[160,95,95,65,65,60];
x = [0:10:n]; y = [0:10:m]; theta = pi/2; s = 1;
lx = length(x); ly = length(y);
min_e = inf;
for i = 1:lx
    for j = 1:ly
        w = [x(i),y(j),theta,s];
        cost = distanceTransform(z,z_part,w,w_part,seq,lF,joint1);
        if cost < min_e
            min_e = cost;
            best_i = i; best_j = j;
        end
    end
end
x = x(best_i); y = y(best_j);
theta = [-pi/2:pi/50:pi/2]; s = [0.2:0.2:0.8 1:0.2:3];
lt = length(theta); ls = length(s);
min_e = inf;
for p = 1:lt
    for q = 1:ls
        w = [x,y,theta(p),s(q)];
        cost = distanceTransform(z,z_part,w,w_part,seq,lF,joint1);
        if cost < min_e
            min_e = cost;
            best_p = p; best_q = q;
        end
    end
end

w = [x,y,theta(best_p),s(best_q)];

x_i = z(1); y_i = z(2); theta_i = z(3); scale_i = z(4);
x_j = w(1); y_j = w(2); theta_j = w(3); scale_j = w(4);
switch w_part
        case 2
            x_ij = x_i-0.5*scale_i*model_len(2)*sin(theta_i) ...
                -0.5*L*cos(theta_i);
            y_ij = y_i+0.5*scale_i*model_len(2)*cos(theta_i) ...
                -0.5*L*sin(theta_i);
            tmp1 = x_j-0.5*scale_j*model_len(2)*sin(theta_j);
            tmp2 = x_j+0.5*scale_j*model_len(2)*sin(theta_j);
            x_ji = bigger(x_ij,tmp1,tmp2);
            tmp1 = y_j+0.5*scale_j*model_len(2)*cos(theta_j);
            tmp2 = y_j-0.5*scale_j*model_len(2)*cos(theta_j);
            y_ji = bigger(y_ij,tmp1,tmp2);
            joint = [x_ji, y_ji];
        case 3
            x_ij = x_i-0.5*scale_i*model_len(3)*sin(theta_i) ...
                +0.5*L*cos(theta_i);
            y_ij = y_i+0.5*scale_i*model_len(2)*cos(theta_i) ...
                +0.5*L*sin(theta_i);
            tmp1 = x_j-0.5*scale_j*model_len(2)*sin(theta_j);
            tmp2 = x_j+0.5*scale_j*model_len(2)*sin(theta_j);
            x_ji = bigger(x_ij,tmp1,tmp2);
            tmp1 = y_j+0.5*scale_j*model_len(2)*cos(theta_j);
            tmp2 = y_j-0.5*scale_j*model_len(2)*cos(theta_j);
            y_ji = bigger(y_ij,tmp1,tmp2);
            joint = [x_ji, y_ji];
    otherwise
end


if w_part == 1|w_part==6
    l = [w(1)-0.5*w(4)*model_len(w_part)*sin(w(3));
    w(2)+0.5*w(4)*model_len(w_part)*cos(w(3));
    w(1)+0.5*w(4)*model_len(w_part)*sin(w(3));
    w(2)-0.5*w(4)*model_len(w_part)*cos(w(3))];
else
    l = [w(1)-0.5*w(4)*model_len(w_part)*sin(w(3)+pi/2);
    w(2)+0.5*w(4)*model_len(w_part)*cos(w(3)+pi/2);
    w(1)+0.5*w(4)*model_len(w_part)*sin(w(3)+pi/2);
    w(2)-0.5*w(4)*model_len(w_part)*cos(w(3)+pi/2)];
end

function a = bigger(b,c,d)
if norm(b-c,1)>norm(b-d,1)
    a = c;
else
    a = d;
end