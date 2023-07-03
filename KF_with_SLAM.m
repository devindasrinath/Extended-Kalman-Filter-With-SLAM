%% EKF SLAM algorithm

%% Step 1 : create basic environment

%map size : 10 x 10
%place land marks in arbitrary points
land_marks_x = [ 1 4 2 -3 -2];
land_marks_y = [ 3 2 -2 3 -2];
figure;
scatter(land_marks_x , land_marks_y , 'r');
xlim([-6,6]);
ylim([-6,6]);
hold on;

%initilalise positions
x_p =0;
y_p=0;
x_r =0;
y_r=0;
u = zeros(12,1);
sigma = [[zeros(2,2)  zeros(2,10)];
         [zeros(10,2) 10000*eye(10,10)]];

%moving velocties and angles
vel = [ 2 2 3 1 1 3 5 1 3 3 ];
theta = [ -30 60 90 120 150 200 240 270 310 30]*(pi/180);

%create noises
noise_var_x = 0.001;
noise_var_y = 0.001;
process_noise_cov = [noise_var_x 0;
                     0 noise_var_y];
sensor_noise_1 = 0.0001;
sensor_noise_2 = 0.0001;
sensor_noises_cov = [sensor_noise_1 0;
                     0 sensor_noise_2];


%% step 2 : create moving model and algorithm
for i = 1:size(vel')
    %planed trajectory
    pre_x_p = x_p;
    pre_y_p = y_p;
    x_p  = pre_x_p + vel(i)*cos(theta(i));
    y_p  = pre_y_p + vel(i)*sin(theta(i));
    plot([pre_x_p x_p] , [pre_y_p y_p],'black');
    
    %real trajectory
    pre_x_r = x_r;
    pre_y_r = y_r;
    x_r  = pre_x_r + vel(i)*cos(theta(i))+ wgn(1,1,noise_var_x,'linear');
    y_r  = pre_y_r + vel(i)*sin(theta(i))+ wgn(1,1,noise_var_y,'linear');
    plot([pre_x_r x_r] , [pre_y_r y_r],'r');
    
    %KF implmentation
    
    %Prediction stage
    F_x = [1 0 0 0 0 0 0 0 0 0 0 0;
           0 1 0 0 0 0 0 0 0 0 0 0];
      
    u_bar = u + F_x'*[vel(i)*cos(theta(i));
                          vel(i)*sin(theta(i))];
   
    G = eye(12);
    
    sigma_bar = G*sigma*G' + F_x'*process_noise_cov*F_x;
    
    %corection step
    %checkfor each land mark
    for j = 0 : 4
        x1 = (land_marks_x(j+1) - x_r);
        y1 = (land_marks_y(j+1) -y_r);
        r = sqrt(x1*x1 + y1*y1);
        beta = atan(y1 / x1);
        if (x1<0) && (y1>0)
            beta = pi - abs(beta);
        elseif (x1<0) && (y1<0)
            beta = pi + abs(beta);
        end
        
        r = r + wgn(1,1,sensor_noise_1,'linear');
        beta = beta + wgn(1,1,sensor_noise_2,'linear');
        
        z=[r;
           beta];
        
        if r >3
            continue;
        end

        %if land mark never seen before do this
        if (u_bar(2*j+3) == 0) && (u_bar(2*j+4) == 0)
            u_bar(2*j+3) = u_bar(1) + r*cos(beta);
            u_bar(2*j+4) = u_bar(2) + r*sin(beta);
        end
        
        delta_x = -u_bar(1) + u_bar(2*j+3);
        delta_y = -u_bar(2) + u_bar(2*j+4);
        
        delta = [delta_x;
                 delta_y];
             
        q = delta'*delta;
        
        beta_1 = atan(delta_y / delta_x);
        if (delta_x<0) && (delta_y>0)
            beta_1 = pi - abs(beta_1);
        elseif (delta_x<0) && (delta_y<0)
            beta_1 = pi + abs(beta_1);
        end
        
        z_i = [ sqrt(q);
              beta_1];
          
        F_x_j = [[1 0;0 1;0 0; 0 0] zeros(4 , 2*j) [0 0;0 0;1 0;0 1] zeros(4, 10 - 2*(j+1))];
        
        H = (1/q)*[(-sqrt(q)*delta_x) (-sqrt(q)*delta_y) (sqrt(q)*delta_x) (sqrt(q)*delta_y); 
                   (delta_y) (-delta_x) (-delta_y) (delta_x)]*F_x_j;
        
        K = sigma_bar*H'*inv(H*sigma_bar*H' + sensor_noises_cov);
        
         u_bar = u_bar + K*(z-z_i);
         sigma_bar = (eye(12,12) - K*H)*sigma_bar;
    end
    u = u_bar;
    sigma = sigma_bar;
   
    %for object circles display
    viscircles([u(1) u(2)] ,100*sqrt(sigma(1,1)*sigma(1,1)+sigma(2,2)*sigma(2,2)));
    
    %for landmarks circle display   
    for k = 1:5
        if(u(2*k+1)~=0 || u(2*k+1)~=0)   
            var_circles(k)  = viscircles([u(2*k+1) u(2*k+2)] , 300*sqrt(sigma(2*k+1,2*k+1)*sigma(2*k+1,2*k+1)+sigma(2*k+2,2*k+2)*sigma(2*k+2,2*k+2)) , 'Color','b');
        end
    end        
    
    pause(1);
 
% uncomment this if you dont need past circles
%     if i <10
%         for k = 1:5
%             if(u(2*k+1)~=0 || u(2*k+1)~=0)
%                 delete(var_circles(k));
%             end
%         end    
%     end
  
end






