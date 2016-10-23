% Simulating lightning 2D case
m = 100;
n = 50;
eta = 5;               % Selectivity
xgrid = 0:m ;           % Initialize domains in x direction 
ygrid = 0:n ;           
xlen = length(xgrid);
ylen = length(ygrid);

A = zeros(xlen*ylen);
b = zeros(xlen*ylen,1);

%for dart = 1:3          % For generating dart leaders
%% We will introduce the set of constraints in the FOR loop itself
for i = 1:xlen
    for j = 1:ylen
        c = (i-1)*ylen + j;                   % Row corresponding to (x,y)  
        c1 = c - ylen;                        % "                  " (x-dx,y) 
        c2 = c + ylen;                        % "                  " (x+dx,y)  
     
        if j == 1                             % Simulate ground here
             A(c,c) = 1;     
             b(c) = 1;
            
        elseif i == 1                         % No flux leakage from left end of system                   
             A(c,[c,c2]) = [1,-1];
             b(c) = 0;
            
        elseif i == xlen                      % No flux leakage from right end of system
             A(c,[c,c1]) = [1,-1];
             b(c) = 0;
             
        elseif j == ylen                      % No flux leakage from top of the system i.e. flux lines don't go back into the clouds
             A(c,[c-1,c]) = [-1, 1];
             b(c) = 0;
             
        else %Condition 2 j's here    
            A(c,[c1, c-1:c+1, c2]) = [1, 1, -4, 1, 1];    % Write the Laplacian here
    
        end
    end
end

%% Specify initial conditions

A((xlen+1)/2*ylen, (xlen+1)/2*ylen - 1) = 0;    % 
A((xlen+1)/2*ylen - 1, :) = 0;                  % Change afterwards
A((xlen+1)/2*ylen - 2, :) = 0;
%A((xlen-1)/2*ylen - 3, :) = 0;

A((xlen+1)/2*ylen - 1, (xlen+1)/2*ylen - 1) = 1;
A((xlen+1)/2*ylen - 2, (xlen+1)/2*ylen - 2) = 1;
%A((xlen-1)/2*ylen - 3, (xlen-1)/2*ylen - 3) = 1;

%% Solve the matrix here                
x = A\b;
B = reshape(x,ylen,[]);   % Reshape vector into nXm matrix. Each column corresponds to unique x
% 
weight = [50 51; 50 50; 50 49; 50 48; 51 48; 52 48; 52 49; 52 50; 52 51;];  % Initialize adjacent domains 

%% Start iterating
while 1
    
    for i = 1:length(weight)
        weight(i,3) = B(weight(i,2), weight(i,1))^eta;    % Store (Potential)^n in third column of w                                         
    end

    weight(:,4) = cumsum(weight(:,3));                    % Store cumulative sum of potentials 
    weight(:,4) = weight(:,4)/weight(length(weight),4);   % Divide by sum of potentials

    random = rand();                                      % Generate random number for determining next step
 
    for i = 1:length(weight)                         % Check whether value lies between two consecutive weights 
        if random < weight(i,4)
            xin = weight(i,1);                       % Store x and y indices  
            yin = weight(i,2);
            A((xin - 1)*ylen + yin, :) = 0;          % Rewrite linear expression for expressing Dirichlet boundary   
            A((xin - 1)*ylen + yin, (xin - 1)*ylen + yin) = 1;
            b((xin - 1)*ylen + yin) = 0;
            break
        end
    end
    
%     sub_weight = weight(weight(:,4) > random, :);
%     sub_weight = sub_weight(1,:);
%      xin = weight(i,1);
%      yin = weight(i,2);
%      A((xin - 1)*ylen + yin, :) = 0;
%      A((xin - 1)*ylen + yin, (xin - 1)*ylen + yin) = 1;
%     
            
    x = A\b ;       
    B = reshape(x,ylen,[]) ;
    B(B < 1e-10 & B > -1e-10) = 0 ;    % Takes care of irregualrities
    %Take 2 readings: one with comment and one sans
    
    %% Updating relevant matrices and book keeping

    if yin == 1 
        break
    else
        ymax = min(yin+1,ylen);         % Checks whether (yin+1) exceeds ylen or not
        xmin = max(xin-1,1);            % "            " (xin-1) is smaller than 1 or not 
        xmax = min(xin+1,xlen);         % "            " (xin+1) exceeds xlen or not
        
        % Calculate indices of domains around the current domain which is
        % being set to zero
        inx = sub2ind(size(B),[yin-1, ymax],[xin, xin]); 
        inx = cat(2,inx,sub2ind(size(B),[yin-1, yin, ymax],[xmax, xmax, xmax]));
        inx = cat(2,inx,sub2ind(size(B),[yin-1, yin, ymax],[xmin, xmin, xmin]));
        
        inx = unique(inx);              % unique(inx, 'rows') for gathering unique rows. Here only unique elements
        inx = inx(B(inx)~=0);           % Only take out those indices which aren't at zero potential
        [y1,x1] = ind2sub(size(B),inx); % Convert back to their original state
        
        add_weight = cat(2,x1',y1',zeros(length(x1),2));    % Add the required new domains 
        weight(i,:) = [];                                   % Erase the domain which was just assigned to 0
        weight = cat(1,weight,add_weight);                  % Append new domains to weight
        
        %mesh(xgrid,ygrid,B);
        %view(0,90);
        %drawnow;
        end
end

%% Plot details
B(B==0) = 5;
B(B~=5) = 0;
mesh(xgrid,ygrid,B);
view(0,90);

title('Path of Lightning Bolt');
xlabel('Ground');
ylabel('Height');
zlabel('Potential');
hold on;
% 
% %% Dart Leaders - Comment out if not needed
% A = zeros(xlen*ylen);               % Add constant charges to produce the dart leaders
% b = reshape(B,[],1);
% %zindex = find(B);                  Finds non zero indices
% b(b==5) = -0.1*dart;                % Here h = 1m. Hence multiplication by h^2 causes no change 
% end