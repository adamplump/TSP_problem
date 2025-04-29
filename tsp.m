% ------------------------------------------------------------
%   INITIAL INFORMATION
% ------------------------------------------------------------

%   Some abbreviations used in variable naming:
%   mpk - cities in order
%   opk - distances in order
%   dn - best found path

%   Program parameter selection by the user

n = 2 + menu('Choose the number of cities to randomly select', '3', '4','5','6','7','8','9','10','11','12','13','14','15');
settings = ["Yes" "No"];
greedy = menu('Do you want to search using the greedy algorithm?',settings);
climbing = menu('Do you want to search using the climbing algorithm?',settings);
depth = menu('Do you want to search using the depth-first algorithm?',settings);

%   Creating a matrix of all distances
m = [65 35; 50 20; 5 65; 95 55; 60 75; 10 30; 80 45; 85 10; 30 70; 30 10; 55 60; 40 40; 80 85; 45 90; 25 45];   %   współrzędne miast
dist_all = zeros(15,15);
for k=1:15                                      
    for w=1:15
        x1 = m(k,1);
        y1 = m(k,2);
        x2 = m(w,1);
        y2 = m(w,2);
        dist_all(w,k)=round(sqrt(((x2-x1)^2)+((y2-y1)^2)));
        if k==w
            dist_all(w,k)=inf;
        end
    end
end

%   Randomly selecting the numbers of chosen cities and creating a distance matrix for them
chosen_city = [1];
chosen_city(2:n) = randperm(14,n-1) + 1;
chosen_dist = zeros(n,n);
for k=1:n
    for w=1:n
        chosen_dist(w,k) = dist_all(chosen_city(w),chosen_city(k));
    end
end

%   Creating a matrix of coordinates for the selected cities
chosen_coord = zeros(n,2);
for k=1:2
    for w=1:n
      chosen_coord(w,k)=m(chosen_city(w),k);
    end
end

% ------------------------------------------------------------
%   GREEDY ALGORITHM
% ------------------------------------------------------------

%   Creating the map
if greedy == 1 
    map('Greedy algorithm', n, chosen_coord);
end

L = chosen_dist;
curr_col = 1;                                     %   Storing the currently checked column
m_city_in_row = zeros(1,n);
m_dist_in_row = zeros(1,n);
for k = 1:n
    A=L(:,curr_col);                              %   A - column being checked in this step   
   [M,I] = min(A);                                %   M - value of the smallest distance, I - index of the smallest distance
   L(curr_col,:) = inf;                           %   Marking the city as already checked
   m_city_in_row(1,k) = I;                        %   Number of the city found in this step
   if greedy == 1
       coord_1 = chosen_coord(curr_col,:);
       coord_2 = chosen_coord(I,:);
       coord_d = coord_2-coord_1;
       quiver(coord_1(1),coord_1(2),coord_d(1),coord_d(2),'LineWidth',1.5,'Autoscale',0,'Color','g')
       axis([0  100    0  100])
   end
   m_dist_in_row(1,k) = chosen_dist(I,curr_col);
   curr_col = I;
end
path_greedy = sum(m_dist_in_row);
if greedy == 1
    disp(['The shortest path found by the greedy algorithm is ', num2str(path_greedy)]);
end

% ------------------------------------------------------------
%   CLIMBING ALGORITHM
% ------------------------------------------------------------

if climbing == 1
    counter = 0;                                 %   Counter for the checked route versions
    counter_best = 0;                            %   Counter for the found shorter routes
    climb_dn = path_greedy;
    climb_mpk = m_city_in_row;
    climb_mpk_temp = climb_mpk;
    climb_mpk_temp2 = climb_mpk;
    step_succ = 1;
    found_better = 0;
    while step_succ == 1 
        for k = 1:n-1
            for w = 2:n-1
                if w>k
                    climb_mpk_temp2 = climb_mpk;
                    climb_mpk_temp2([k w]) = climb_mpk_temp2([w k]);
                    for i = 1:n
                        if i==1
                            climb_opk(1,i) = chosen_dist(1,climb_mpk_temp2(i));
                        else
                            climb_opk(1,i) = chosen_dist(climb_mpk_temp2(i-1),climb_mpk_temp2(i));
                        end
                    end
                    climb_dc = sum(climb_opk);
                    if climb_dc < climb_dn 
                        climb_dn = climb_dc;
                        climb_mpk_temp = climb_mpk_temp2;
                        found_better = 1;
                    end
                    counter = counter+1;
                end
            end   
        end
        if found_better == 1
            found_better = 0;
            counter_best = counter_best+1;
            step_succ = 1;
            map(['Climbing algorithm - step ', num2str(counter_best)], n, chosen_coord);
            for i = 1:n
                if i==1
                    coord_1 = chosen_coord(1,:);
                    coord_2 = chosen_coord(climb_mpk(i),:);
                    coord_d = coord_2 - coord_1;
                    quiver(coord_1(1),coord_1(2),coord_d(1),coord_d(2),'LineWidth',4,'Autoscale',0,'Color','c')
                else
                    coord_1 = chosen_coord(climb_mpk(i-1),:);
                    coord_2 = chosen_coord(climb_mpk(i),:);
                    coord_d = coord_2 - coord_1;
                    quiver(coord_1(1),coord_1(2),coord_d(1),coord_d(2),'LineWidth',4,'Autoscale',0,'Color','c')
                end
            end
            for i = 1:n
                if i==1
                    coord_1 = chosen_coord(1,:);
                    coord_2 = chosen_coord(climb_mpk_temp(i),:);
                    coord_d = coord_2 - coord_1;
                    quiver(coord_1(1),coord_1(2),coord_d(1),coord_d(2),'LineWidth',1.5,'Autoscale',0,'Color','k')
                else
                    coord_1 = chosen_coord(climb_mpk_temp(i-1),:);
                    coord_2 = chosen_coord(climb_mpk_temp(i),:);
                    coord_d = coord_2 - coord_1;
                    quiver(coord_1(1),coord_1(2),coord_d(1),coord_d(2),'LineWidth',1.5,'Autoscale',0,'Color','k')
                end
            end
            climb_mpk = climb_mpk_temp;
        else
            step_succ = 0;
        end      
    end
    if counter_best == 0
        disp(['The initial path for the climbing algorithm was ', num2str(path_greedy),' and it did not find a shorter path']);
    elseif counter_best == 1
        disp(['The climbing algorithm performed ', num2str(counter_best),' step, and the shortest found path is ', num2str(climb_dn)]);
    elseif counter_best >= 5
        disp(['The climbing algorithm performed ', num2str(counter_best),' steps, and the shortest found path is ', num2str(climb_dn)]);
    else
        disp(['The climbing algorithm performed ', num2str(counter_best),' steps, and the shortest found path is ', num2str(climb_dn)]);
    end
end

% ------------------------------------------------------------
%   DEPTH-FIRST ALGORITHM
% ------------------------------------------------------------

if depth == 1
    range = chosen_city(2:n);
    for k = 1:n
        depth_start_order(1,k) = chosen_city(m_city_in_row(k));
    end
    
    [ich,depth_best_way,depth_mpk] = depth_alg(depth_start_order,dist_all,path_greedy,range);

%   Drawing maps  
    map('Depth-first algorithm (brute force)', n, chosen_coord);
    for i = 1:n
                if i==1
                    coord_1 = m(1,:);
                    coord_2 = m(depth_mpk(i),:);
                    coord_d = coord_2-coord_1;
                    quiver(coord_1(1),coord_1(2),coord_d(1),coord_d(2),'LineWidth',1.5,'Autoscale',0,'Color','b')
                else
                    coord_1 = m(depth_mpk(i-1),:);
                    coord_2 = m(depth_mpk(i),:);
                    coord_d = coord_2-coord_1;
                    quiver(coord_1(1),coord_1(2),coord_d(1),coord_d(2),'LineWidth',1.5,'Autoscale',0,'Color','b')
                end
    end
    disp(['The shortest path found by the depth-first algorithm (brute force) is ', num2str(depth_best_way)]);
end

% ------------------------------------------------------------
%   FUNCTION DEFINITIONS
% ------------------------------------------------------------

function map(name_map, n, chosen_coord)
figure('name',name_map);
for w=1:n
    plot(chosen_coord(w,1),chosen_coord(w,2),'ok','LineWidth',5);     % Marking the selected cities on the plot
    if w==1
        plot(chosen_coord(w,1),chosen_coord(w,2),'or','LineWidth',7);   % Highlighting the starting point
    end
    hold on
end
axis([0  100    0  100])
end

function [odp,depth_dn,depth_mpk] = depth_alg(depth_start_order,dist_all,path_greedy,range)
persistent curr_route
persistent brute_dn
persistent depth_mpk_temp
persistent quant_cities
if isempty(brute_dn)
    brute_dn = path_greedy;
end
n = numel(range);
if isempty(quant_cities)
    quant_cities = n+1;
end
if isempty(depth_mpk_temp)
    depth_mpk_temp = depth_start_order;
end
if n <= 1
    odp = range;
    curr_route(n) = range;
    for k = 1:(numel(curr_route)+1)
        if k == 1
            curr_path(1,k) = dist_all(curr_route(k),1);
        elseif k == numel(curr_route)+1
            curr_path(1,k) = dist_all(1,curr_route(k-1));
        else
            curr_path(1,k) = dist_all(curr_route(k),curr_route(k-1));
        end
    end
    if sum(curr_path) < brute_dn
        brute_dn = sum(curr_path);
        depth_mpk_temp = curr_route;
        depth_mpk_temp(quant_cities) = 1;
    end
else
    odp = zeros(factorial(n),n);
    f = factorial(n-1);
    L = 1:f;
    for i = n:-1:1
        odp(L,1) = range(i);
        curr_route(n) = range(i);
        odp(L,2:n) = depth_alg(depth_start_order,dist_all,path_greedy,range(setdiff(1:n,i)));
        L = L + f;
    end
end
depth_dn = brute_dn;
depth_mpk = depth_mpk_temp;
end
