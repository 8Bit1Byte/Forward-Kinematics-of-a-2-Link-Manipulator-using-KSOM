%% POINTS TRACKING

u1 = [ 0 1.414; 1.414 0; 1 1; -1 -1; 0.2 0.8]';

v = zeros(5, 4);
for m = 1:size(u1, 2)
    u = u1(:, m);
    for j = 1:100
        dist(j) = norm(u-w_g(:, :, j));
    end
    
    [~, win_val] = min(dist);
    win = [lx(win_val), ly(win_val)];
    d = repmat(win', 1, 100)-lattice;
    H_g = exp(-(sum(d.^2))/(2*(sig_f^2)));

    % Corse Action
    s = sum(H_g); s2 = 0; s3 = 0;
    
    for k = 1:100
        s1 = H_g(k)*(th_g(:,:,k)+A_g(:,:,k)*(u-w_g(:,:,k)));
        s2 = s2 + s1;
    end      
    
    theta = s2/s;
    x = l1*cos(theta(1)) + l2*cos(theta(1) + theta(2));
    y = l1*sin(theta(1)) + l2*sin(theta(1) + theta(2));
    % Tracked Point
    v(m, 1) = x; v(m, 2) = y; 
    v(m, 3) = theta(1)*180/pi; v(m, 4) = theta(2)*180/pi;
end
array2table(v,...
    'VariableNames', {'X', 'Y', 'theta_1', 'theta_2'},...
    'RowNames', ...
    {'(0.000,1.414)','(1.414,0.000)',...
    '(1.000,1.000)','(-1.00,-1.00)','(0.200,0.800)'})