%% KSOM TRANING NETROWK
% -----------------------

clc; clear; close all;

% Inilialization of Model Params
sig_i= 2.5; sig_f = 0.01;
etaw_i = 1; etaw_f = 0.05;
etaA_i = 0.9; etaA_f = 0.9;

% Link lenghts
l1 = 1; l2 = 1;

A_g = 0.1*rand(2, 2, 100);
w_g = 0.1*rand(2, 1, 100);
th_g = 0.1*rand(2, 1, 100);

% 2d lattice formation of size 10x10
[lx, ly] = ind2sub([10, 10], 1:100);
lattice = [lx; ly]; iterations = 6000;

% Iterations and update
for i = 1:iterations
    th1 = (rand - 0.5)*2*pi; th2 = (rand - 0.5)*2*pi;
    x = l1*cos(th1) + l2*cos(th2 + th1);
    y = l1*sin(th1) + l2*sin(th2 + th1);

    u = [x; y];
    for j = 1:100
        dist(j) = norm(u-w_g(:,:,j));
    end

    [~, win_val] = min(dist);

    % Winning Neuron
    win = [lx(win_val), ly(win_val)];
    sig(i) = sig_i*((sig_f/sig_i)^(i/iterations));
    eta_wg(i) = etaw_i*((etaw_f/etaw_i)^(i/iterations));
    eta_Ag(i) = etaA_i*((etaA_f/etaA_i)^(i/iterations));
    d = repmat(win', 1, 100)-lattice;
    H_g = exp( - (sum(d.^2))/(2*(sig(i)^2)));

    % Coarse action
    s = sum(H_g); s2 = 0; s3 = 0;
    for k = 1:100
        s1 = H_g(k)*(th_g(:,:,k)+A_g(:,:,k)*(u-w_g(:,:,k)));
        s2 = s2 + s1;
    end

    th_o = s2/s;
    x_o = l1*cos(th_o(1)) + l2*cos(th_o(1) + th_o(2));
    y_o = l1*sin(th_o(1)) + l2*sin(th_o(1) + th_o(2));
    v_o = [x_o; y_o];

    % Fine action
    for k = 1:100
        s4 = H_g(k)*(A_g(:,:,k)*(u-v_o));
        s3 = s3 + s4;
    end

    th_1 = th_o + s3/s;
    x_1 = l1*cos(th_1(1)) + l2*cos(th_1(1) + th_1(2));
    y_1 = l1*sin(th_1(1)) + l2*sin(th_1(1) + th_1(2));
    v_1 = [x_1; y_1];

    % Update equtions
    del_v = v_1-v_o;
    del_th = th_1-th_o;
    s5 = 0; s7 = 0;

    for k = 1:100
        s6 = H_g(k)*(th_g(:,:,k)+A_g(:,:,k)*(v_o-w_g(:,:,k)));
        s5 = s5 + s6;
    end

    for t = 1:100
        deltheta_g(:,:,t) = (H_g(t)/s)*(th_o-(s5/s));
    end

    for k = 1:100
        s8 = H_g(k)*(A_g(:,:,k)*del_v);
        s7 = s8 + s7;
    end

    for t = 1:100
        deltaA_g(:,:,t) = (H_g(t)/(s*norm(del_v)^2))*(del_th-s7/s).*(del_v');
        w_g(:,:,t) = w_g(:,:,t) + eta_wg(i)*H_g(t)*(u-w_g(:,:,t));
        th_g(:,:,t) = th_g(:,:,t) + eta_Ag(i)*deltheta_g(:,:,t);
        A_g(:,:,t) = A_g(:,:,t) + eta_Ag(i)*deltaA_g(:,:,t);
    end
end

% Plot final Weights
figure(1); hold on;
for t = 1:100
    plot(w_g(1,1,t), w_g(2,1,t), '*');
end
