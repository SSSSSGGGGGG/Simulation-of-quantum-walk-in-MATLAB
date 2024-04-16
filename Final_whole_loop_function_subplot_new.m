clc
close all
clear all

Input_pol=[1;0];

phi=170; % phase difference between fast and slow axises
phi_=(phi/180)*pi;
HWP=[1,0;0,exp(-i*phi_)]; 
 
a=25;  % angle of HWP, change it! 
theta=(a/180)*pi;% rotated angle
c_theta=cos(theta); %cos of theta
s_theta=sin(theta); %sin of theta
R=[c_theta,s_theta;-s_theta,c_theta]; %rotate matrix 
R_nag=[c_theta,-s_theta;s_theta,c_theta];
HWP_R=R_nag*HWP*R; % HWP with rotation of theta

n=5; % quantum walker step

if n==1 % if the n=1 
   W_output = calculate(n,a); % use the result of function when n=1 
   i=1;       % parameters to calculate the indensity percentage after walking.
   Total_i=0; % the total intensity
   while i<=2   % row
       j=1;         % column
       while j<=2
         Total_i= Total_i+(W_output(j,i))^2;
       % so far we have the total intensity as denominator
        j=j+1;
       end       
       i=i+1;       
   end
   i=1; % initial i to calculate sub-intensity act as column
   I_1=0;   % the intensity of each vector (one column)
   while i<=2 % column
       k=1; % row
       sum_2=0;
       sum_1=0; 
       if k==1 && k<=2
        sum_1=sum_1+(W_output(k,i))^2;
        k=k+1;
        if k==2
        sum_2=sum_2+(W_output(k,i))^2;
        end
        I_1(i)=(sum_1+sum_2)/Total_i;  % percentage of intensity
       end
       k=k+1;
       i=i+1; 
   end
   W_output;
    % figure('Position',[500,450,400,400]);
    % figure
    subplot(1,2,1);
    polellip(W_output(1:2,1)); title("intensity="+I_1(1)*100+' %')
    % figure.position(1:4)=[1,500,300,300];
    xlim ([-1 1])
    ylim ([-1 1])
    % figure('Position',[500,50,400,400]);
    subplot(1,2,2);
    polellip(W_output(1:2,2)); title("intensity="+I_1(2)*100+' %')
    xlim ([-1 1])
    ylim ([-1 1])
    
   % [tau_0,epsilon_0,ar_0,rs_0] = polellip(W_output);
   
else         % when n is larger than 1
    W_output_2 = calculate(1,a);
    C_out_pol = HWP_R*W_output_2
    while 1<n 
    [x,y]=size(C_out_pol);
    W_output_2=size(x,y+1);
    W_output_2=C_out_pol;
    j=y;
        while 0<j    % judge the value along y_axis, if it is not 0, and move it to the next position
             if abs(C_out_pol(2,j))>0
                W_output_2(2,j+1)=C_out_pol(2,j);
                W_output_2(2,1)=0; % this postion is always 0, due to there is no preceeded position.
             else
                W_output_2(2,j)=0; % to check the error
             end
           
            j=j-1;

        end
        
    W_output_2
    C_out_pol = HWP_R*W_output_2
    
    n=n-1;
    end
    [x_1,y_1]=size(W_output_2);
    i=1;       % parameters to calculate the indensity percentage after walking.
   Total_i=0; % the total intensity
   while i<=x_1   % row
       j=1;         % column
       while j<=y_1
         Total_i= Total_i+(W_output_2(i,j))^2;
       % so far we have the total intensity as denominator
        j=j+1;
       end       
       i=i+1;       
   end
    i=1; % initial i to calculate sub-intensity act as column
   I_1=0;   % the intensity of each vector (one column)
   while i<=y_1 % column
       k=1; % row
       sum_2=0;
       sum_1=0; 
       if k==1 && k<=2
        sum_1=sum_1+(W_output_2(k,i))^2;
        k=k+1;
        if k==2
        sum_2=sum_2+(W_output_2(k,i))^2;
        end
        I_1(i)=(sum_1+sum_2)/Total_i;  % percentage of intensity
       end
       k=k+1;
       i=i+1; 
   end
    
    syms r  % number of row
    eqn=r*r==y_1; % y_1 means how many figures you want to show
    S=solve(eqn, r); 
    z=round(double(S(sign(S)==1))) % to calculate how many rows we should have
            for m=1:y_1 
            [tau(m),epsilon(m),ar(m),rs(m)]=polellip(W_output_2(1:2,m));
            subplot(1,y_1,m);polellip(W_output_2(1:2,m));title(rs(m)+ " I="+real(I_1(m)*100)+' %')

            xlim ([-1 1])
            ylim ([-1 1])

            end
    
end






function W_output = calculate(n,a)
    Input_pol=[1;0]
    HWP=[1,0;0,-1]; 
    % a=22.5;  % angle of HWP, change it!
    theta=(a/180)*pi;% rotated angle, later it will be changed
    c_theta=cos(theta); %cos of theta
    s_theta=sin(theta); %sin of theta
    R=[c_theta,s_theta;-s_theta,c_theta]; %rotate matrix 
    R_nag=[c_theta,-s_theta;s_theta,c_theta];
    HWP_R=R_nag*HWP*R; % HWP with rotation of theta
    % n=1;
    C_out_pol=HWP_R*Input_pol % after coin
    W_output=size(2,2);
    W_output=C_out_pol;
    while 0<n 
        if abs(C_out_pol(2,n))>0
                W_output(2,n+1)=C_out_pol(2,n);
                W_output(2,1)=0;
        else
           W_output(2,n+1)=0; 
        end
        n=n-1;
    end
    W_output
end
