function D_f = distanceTransform(z,z_part,w,w_part,seq,lF,joint)
% z is known and w is unknown, distance transform is a minimization
% problem to get optimal w

%% parameters
w_theta = 50; w_scale = 10; theta_ij = 0; scale_ij = 1; W_ij = diag([3,3]); L = 130;
% w_theta = 0; w_scale = 0; theta_ij = 0; scale_ij = 1; W_ij = diag([3,3]); L = 130;


%% Transformation of variables
% [x_i, y_i, theta_i, scale_i] = z([1 2 3 4]);
% [x_j, y_j, theta_j, scale_j] = w([1 2 3 4]);
x_i = z(1); y_i = z(2); theta_i = z(3); scale_i = z(4);
x_j = w(1); y_j = w(2); theta_j = w(3); scale_j = w(4);
% joint points [x_ij, y_ij]  [x_ji, y_ji]
% for w_part, 2 = left upper arm, 3 = right upper arm, 6 = head
model_len=[160,95,95,65,65,60];
if z_part == 1
    switch w_part
        case 2
            x_ij = x_i-0.5*scale_i*model_len(2)*sin(theta_i) ...
                -0.5*L*cos(theta_i);
            y_ij = y_i+0.5*scale_i*model_len(2)*cos(theta_i) ...
                -0.5*L*sin(theta_i);
            tmp1 = x_j-0.5*scale_j*model_len(2)*sin(theta_j);
            tmp2 = x_j+0.5*scale_j*model_len(2)*sin(theta_j);
            x_ji = smaller(x_ij,tmp1,tmp2);
            tmp1 = y_j+0.5*scale_j*model_len(2)*cos(theta_j);
            tmp2 = y_j-0.5*scale_j*model_len(2)*cos(theta_j);
            y_ji = smaller(y_ij,tmp1,tmp2);
        case 3
            x_ij = x_i-0.5*scale_i*model_len(3)*sin(theta_i) ...
                +0.5*L*cos(theta_i);
            y_ij = y_i+0.5*scale_i*model_len(2)*cos(theta_i) ...
                +0.5*L*sin(theta_i);
            tmp1 = x_j-0.5*scale_j*model_len(2)*sin(theta_j);
            tmp2 = x_j+0.5*scale_j*model_len(2)*sin(theta_j);
            x_ji = smaller(x_ij,tmp1,tmp2);
            tmp1 = y_j+0.5*scale_j*model_len(2)*cos(theta_j);
            tmp2 = y_j-0.5*scale_j*model_len(2)*cos(theta_j);
            y_ji = smaller(y_ij,tmp1,tmp2);
        case 6
            x_ij = x_i-0.5*scale_i*model_len(6)*sin(theta_i);
            y_ij = y_i+0.5*scale_i*model_len(6)*cos(theta_i);
            x_ji = x_j+0.5*scale_j*model_len(6)*sin(theta_j);
            y_ji = y_j-0.5*scale_j*model_len(6)*cos(theta_j);
    end
else
    x_ij = joint(1); y_ij = joint(2);
    tmp1 = x_j-0.5*scale_j*model_len(2)*sin(theta_j);
    tmp2 = x_j+0.5*scale_j*model_len(2)*sin(theta_j);
    x_ji = smaller(x_ij,tmp1,tmp2);
    tmp1 = y_j+0.5*scale_j*model_len(2)*cos(theta_j);
    tmp2 = y_j-0.5*scale_j*model_len(2)*cos(theta_j);
    y_ji = smaller(y_ij,tmp1,tmp2);
end


theta_i = w_theta*(theta_i-theta_ij/2);
theta_j = w_theta*(theta_j-theta_ij/2);
scale_i = w_scale*(log(scale_i)-log(scale_ij)/2);
scale_j = w_scale*(log(scale_j)-log(scale_ij)/2);
xy_i = W_ij*([x_ij, y_ij]'); 
xy_j = W_ij*([x_ji, y_ji]'); 
T_ij = [theta_i, scale_i, xy_i'];
T_ji = [theta_j, scale_j, xy_j'];

%% cost function
f_w = match_energy_cost(w,w_part,seq,lF);

D_f = norm(T_ij-T_ji,1)+2*f_w;

function a = smaller(b,c,d)
if norm(b-c,1)<norm(b-d,1)
    a = c;
else
    a = d;
end




