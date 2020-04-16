clc; t = []; syms t;

% PATH PARAMS
len = 100; % number of data points to generate velocities for along curve
b = 0.34; % speed coefficient
bt_max = 3.2 - 0.15; % endpoint of path (lowered to prevent falling off bridge)
R = 4*[0.396*cos(2.65*(b*t+1.4)); -0.99*sin(b*t+1.4); 0]; % path function in terms of t
filename = 'BODpath.mat'; % filename to save pathdata under

% GENERATE PATHDATA
calculate(R, bt_max/b, len, filename);

function calculate(R, t_max, len, filename)
% Calculates the right & left wheelspeeds for len number of data points to 
% execute path R over time, saves them under filename to be executed later.
    t = []; syms t;
    d = .235; % wheelbase width of neato in meters

    [left, right, That] = calc_functions(R, d);
    pathdata = generate_vel(R, left, right, len, t_max);
    save(filename, 'pathdata', 'R', 'That')
    disp('Done with calculations')

    function [left, right, That] = calc_functions(R, d)
    % Symbolically calculates the functions for the right and left
    % wheelspeeds for path 'R' with wheelbase width 'd'. Also saves the
    % unit tangent vector for getting the initial heading of the neato.
        T=diff(R,t); % linear velocity vector
        That=(T./norm(T));
        dThat=diff(That,t);
        Bhat=(cross(That,dThat));
        w = Bhat; % angular velocity vector

        left = norm(T) - norm(w)*d/2;
        right = norm(T) + norm(w)*d/2;
    end

    function pathdata = generate_vel(R, left, right, len, t_max)
    % Discretizes into data points by calculating 'len' number of evenly 
    % spaced points over time along curve 'R' and wheelspeed functions 
    % 'right' and 'left' from time 0 to 't_max'.
        t_num =linspace(0,t_max,len);
        r_num =zeros(len,3);
        left_num=zeros(len,1);
        right_num=zeros(len,1);

        for n=1:len
            r_num(n,:)=double(subs(R,t,t_num(n)));
            left_num(n)=double(subs(left,t,t_num(n)));
            right_num(n)=double(subs(right,t,t_num(n)));
        end
        pathdata = [transpose(t_num), r_num(:,1:2), left_num, right_num];
    end
end