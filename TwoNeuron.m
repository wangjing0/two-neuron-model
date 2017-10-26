%% Two inhibitory neurons with various level of common input
% for generating the nueral trajectories, speed and energy profiles, as shown in Figure 6c and 6b
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% by Jing Wang 1/29/2017  jingwang.physics(a)gmail.com
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clc ; clear;  close all
%% PARAMETERS
        Ntrial = 1;
         Tau   = 100;
            Dt = 10;
            dt = (Dt/Tau); % in unit of tau
            Nt = 200;
        tictoc = ceil(800/Dt);
                         
        Ntheta = 4;
         Theta = linspace(.6,.8,Ntheta);
   Theta_Noise = 0.001;
     x0_jitter = .05;
             W = 6.0;

%% Initialization
     xmin  = 0;
     xmax  = 1.2;
     x     = linspace(xmin,xmax,1000);
     cc    = flipud(linspecer(Ntheta*10-9));
     Color = cc(1:10:end,:);

         X = nan(Ntheta,2,Nt);
         V = nan(Ntheta,2,Nt);
     Speed = nan(Ntheta,Nt);
    Energy = nan(Ntheta,Nt);
     Traj  = nan(Ntheta,Ntrial,Nt,2);
    
figure('Position',[100 500 1200 300]);
subplot(131)
xlabel('u') ;ylabel('v')
xlim([xmin xmax]); ylim([xmin xmax])
set(gca, 'XTick', [0 1], 'YTick', [0 1]);
drawnow; hold on
subplot(132)
xlabel('u'); ylabel('v') ;zlabel('Speed')
xlim([xmin xmax]); ylim([xmin xmax])
set(gca, 'XTick', [], 'YTick', [],'ZTick', []);
view([-45 10]); box off
drawnow; hold on
subplot(133)
xlabel('u'); ylabel('v') ;zlabel('Energy')
xlim([xmin xmax]); ylim([xmin xmax])
set(gca, 'XTick', [], 'YTick', [],'ZTick', []);
view([-45 10]);
drawnow; hold on

%% trajectory
for ii=1:Ntheta
    
    clear  y xx yy u v   
    y=1./(1+exp(-1.0 * W * (Theta(ii) - x)));
    
    [u0,v0]=findIntersect(x(x>=0.5 & x<0.7 ),y(x>=0.5  & x<0.7),y,x);


    subplot(131)  
        plot(x,y,'Color',Color(ii,:),'LineWidth',0.2); hold on 
        plot(y,x,'--','Color',Color(ii,:),'LineWidth',0.2); hold on
        plot(u0,v0,'o','MarkerEdgeColor',Color(ii,:)); hold on; drawnow 
   subplot(132)
        plot3(u0,v0,0,'o','MarkerEdgeColor',Color(ii,:)); hold on; drawnow
    subplot(133)
        plot3(u0,v0,0,'o','MarkerEdgeColor',Color(ii,:)); hold on; drawnow
    
    [u1,v1]=findIntersect(x,y,x,x + x0_jitter);
    [u2,v2]=findIntersect(y,x,x,x + x0_jitter);
    u0=mean([u1,u2]); v0=mean([v1,v2]);
   
    for j=1:Ntrial
        
        v = v0 ;
        u = u0 ;
        
        for t=1:Nt
            Traj(ii,j,t,:)=[u v];
            [Vx_,~]= findIntersect(linspace(xmin,xmax,1e3),v*ones(1,1e3),y,x);
            [~,Vy_]= findIntersect(u*ones(1,1e3),linspace(xmin,xmax,1e3),x,y);
 
            ksi=randn; % correlated input noise to both u and v
            
            Vx = - u + Vx_ + Theta_Noise*ksi;
            Vy = - v + Vy_ + Theta_Noise*ksi;
            
            v = v + Vy.*dt;
            u = u + Vx.*dt;
            if j==1
                 X(ii,:,t) = [u;v];
                 V(ii,:,t) = [Vx;Vy];
                Speed(ii,t)= norm([Vx,Vy],2);
            end
            
            subplot(131)
            if ~mod(t,5)
            plot(u,v,'.','Color',Color(ii,:)); hold on; drawnow
            end
            
            if t ==  tictoc  % mark the interval
                plot(u,v,'kd'); hold on  ; drawnow
            end           
        end
        plot(u,v,'o','MarkerEdgeColor',Color(ii,:),'MarkerFaceColor',Color(ii,:)); hold on ; drawnow
    end
    
    Energy(ii,1) = 0;
    Energy(ii,2:end) = -1.0 * cumsum(diff((X(ii,1,:))).*((V(ii,1,2:end))) + diff((X(ii,2,:))).*((V(ii,2,2:end)))); 
   
    subplot(132)
    plot3(squeeze(X(ii,1,:)),squeeze(X(ii,2,:)), squeeze(Speed(ii,:)) ,'.','Color',Color(ii,:)); hold on; drawnow
    subplot(133)
    plot3(squeeze(X(ii,1,:)),squeeze(X(ii,2,:)), squeeze(Energy(ii,:)) ,'.','Color',Color(ii,:)); hold on; drawnow
end

