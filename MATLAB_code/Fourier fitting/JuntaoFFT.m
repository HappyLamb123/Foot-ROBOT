function [alphaz_co,alphax_co,mn_az,mn_ax,MN_shift_z,MN_shift_x,filtaz,filtax] = JuntaoFFT(n) 
% function [filtaz,filtax] = JuntaoFFT(n) 
%%%%%%%%% return the alphaz coefficients and alphax coefficients
%%%%%%%%% the argument is the fitting orders
    %Load the data from the mat files\
%     clear
%     clc
%     n = 2
    load data1
    load data2
    load data3
    %Average the three data sets
    jnk=(data_1+data_2+data_3)/3;
   
    
    
    %Extract the two alpha matrices
    jnkax=jnk(10:end,2:end);
    jnkaz=jnk(2:8,2:end);
    
    %Extract the gamma and beta
    gamma = jnk(1,2:end);
    beta = jnk(2:8,1);
    left_range = 2:7;
    right_range = 8:13;
    %Make the matrices periodic 
    % jnkax2=[jnkax, fliplr(jnkax)];
    % jnkaz2=[jnkaz, fliplr(jnkaz)];
%     jnkax2=[-fliplr(jnkax(:,1:7)),jnkax, -fliplr(jnkax(:,8:14))];
%     jnkaz2=[fliplr(jnkaz(:,1:7)), jnkaz, fliplr(jnkaz(:,8:14))];

    jnkax2=[-fliplr(jnkax(:,2:7)),jnkax, -fliplr(jnkax(:,8:13))];
    
    jnkaz2=[fliplr(jnkaz(:,2:7)), jnkaz, fliplr(jnkaz(:,8:13))];
    
    gamma_flip = [fliplr(-pi-gamma(left_range)) gamma fliplr(pi-gamma(right_range))];
    gamma_interval = (gamma_flip(1,end)-gamma_flip(1,1))/25;
    gamma_range = gamma_flip(1,1):gamma_interval:gamma_flip(1,end);
    
%     [X,Y] =meshgrid(-pi:pi*2/25:pi,-pi/2:pi/6:pi/2);
    [X,Y] =meshgrid(gamma_range,pi/2:-pi/6:-pi/2);
    
    
    [Gamma,Beta]=meshgrid(gamma_flip,beta);
    % interpolate the data to make the data evenly distributed in terms of
    % gamma
    x_interp = interp2(Gamma,Beta,jnkax2,X,Y,'linear');
%     x_interp(:,1)=2*x_interp(:,2)-x_interp(:,3);
%     x_interp(:,end)=2*x_interp(:,26-1)-x_interp(:,26-2);
    z_interp = interp2(Gamma,Beta,jnkaz2,X,Y);
%     z_interp(:,1)=2*z_interp(:,2)-z_interp(:,3);
%     z_interp(:,end)=2*z_interp(:,26-1)-z_interp(:,26-2);
    % jnkax2=jnkax;
    % jnkaz2=jnkaz;
    %Plot the raw data
    
    figure(1)
    subplot(2,2,1);
%     imagesc([-pi,pi],[-pi/2,pi/2],jnkaz2)
    imagesc([-pi,pi],[-pi/2,pi/2],z_interp)
    axis equal
    axis([-pi,pi,-pi/2,pi/2]);
    colorbar
    xlabel('$\gamma$')
    ylabel('$\beta$');
    title('$\alpha_z$');
    set(gca,'FontSize',20)
    subplot(2,2,2);
%     imagesc([-pi,pi],[-pi/2,pi/2],jnkax2)
    imagesc([-pi,pi],[-pi/2,pi/2],x_interp)  
    axis equal
    axis([-pi,pi,-pi/2,pi/2]);
    colorbar
    xlabel('$\gamma$')
    ylabel('$\beta$');
    title('$\alpha_x$');
    set(gca,'FontSize',20)
    %Fourier transform the data
%     fftaz=fftshift(fft2(jnkaz2));
%     fftax=fftshift(fft2(jnkax2));
    fftaz=fftshift(fft2(z_interp));
    fftax=fftshift(fft2(x_interp));
    %Make a mask to filter the data 
    mask=zeros(size(fftaz));
    %Include n lowest frequency components
%     n=2;
    mid_gamma = 14; %15
    mask((4-n):(n+4),(mid_gamma-n):(mid_gamma+n))=1;

    % mask((4-n):(n+4),(8-n):(8+n))=1;
    %Apply the filter
%     fft2(jnkaz2)
    fftazfilt=fftaz.*mask;
    fftaxfilt=fftax.*mask;
    %Transform the arrays back to real space
%     filtaz=real(ifft2(ifftshift(fftazfilt)));
%     filtax=real(ifft2(ifftshift(fftaxfilt)));
    filtaz=ifft2(ifftshift(fftazfilt));
    filtax=ifft2(ifftshift(fftaxfilt));
    %Plot the filtered arrays
    subplot(2,2,3);
    imagesc([-pi,pi],[-pi/2,pi/2],filtaz);
    axis equal
    axis([-pi,pi,-pi/2,pi/2]);
    colorbar
    xlabel('$\gamma$')
    ylabel('$\beta$');
    title('$\alpha_z$ filtered');
    set(gca,'FontSize',20)
    subplot(2,2,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot error map%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imagesc([-pi,pi],[-pi/2,pi/2],filtax);
    axis equal
    axis([-pi,pi,-pi/2,pi/2]);
    colorbar
    xlabel('$\gamma$')
    ylabel('$\beta$');
    title('$\alpha_x$ filtered');
    set(gca,'FontSize',20)
    figure
    subplot(1,2,1);
    
    imagesc([-pi,pi],[-pi/2,pi/2],filtax-x_interp);
    axis equal
    axis([-pi,pi,-pi/2,pi/2]);
    colorbar
    xlabel('$\gamma$')
    ylabel('$\beta$');
    title(['$\alpha_x$ error Order:',num2str(n)]);
    set(gca,'FontSize',20)
    caxis([-0.1,0.1]);
    subplot(1,2,2);
    
    imagesc([-pi,pi],[-pi/2,pi/2],filtaz-z_interp);
    axis equal
    axis([-pi,pi,-pi/2,pi/2]);
    colorbar
    xlabel('$\gamma$')
    ylabel('$\beta$');
    title(['$\alpha_z$ error  Order:',num2str(n)]);
    set(gca,'FontSize',20)
    caxis([-0.1,0.1]);
    colormap([[(0:127)/127; (0:127)/127; ones(1,128)],[ones(1,128); (127:-1:0)/127; (127:-1:0)/127]]')
%     alphaz_co=fftazfilt((4-n):(n+4),(15-n):(15+n))/28/7;
%     alphax_co=fftaxfilt((4-n):(n+4),(15-n):(15+n))/28/7;
%     aaz = fft2(jnkaz2)/28/7;
%     aax = fft2(jnkax2)/28/7;
%     aaz = ifftshift(fftazfilt)/28/7;
%     aax = ifftshift(fftaxfilt)/28/7;
    aaz = ifftshift(fftazfilt)/26/7;
    aax = ifftshift(fftaxfilt)/26/7;
    alphaz_co=[];
    alphax_co=[];
    mn_az = [];
    mn_ax = [];
    MN_shift_z= [];
    MN_shift_x= [];
    M_z = 0;
    N_z = 0;
    M_x = 0;
    N_x = 0;
    for i = 1:7
%         for j = 1:28
        for j = 1:26 %28 
            if abs(aaz(i,j))>0.00005
                alphaz_co=[alphaz_co,aaz(i,j)];
                mn_az = [mn_az;i,j];
%                 if i > 4
%                     M_z = i - 8;
%                 else
%                     M_z = i-1;
%                 end
%                 if j > 15
%                     N_z = j - 29;
%                 else 
%                     N_z = j-1;
%                 end
                if i > 4
                    M_z = i - 8;
                else
                    M_z = i-1;
                end
                if j > mid_gamma
                    N_z = j - (mid_gamma*2-1);
                else 
                    N_z = j-1;
                end
                MN_shift_z= [MN_shift_z;M_z,N_z];
            end
            if abs(aax(i,j))>0.00001
                alphax_co=[alphax_co,aax(i,j)];
                mn_ax = [mn_ax;i,j];
%                 if i > 4
%                     M_x = i - 8;
%                 else
%                     M_x = i-1;
%                 end
%                 if j > 15
%                     N_x = j - 29;
%                 else 
%                     N_x = j-1;
%                 end
                if i > 4
                    M_x = i - 8;
                else
                    M_x = i-1;
                end
                if j > mid_gamma
                    N_x = j - (2*mid_gamma-1);
                else 
                    N_x = j-1;
                end
                MN_shift_x= [MN_shift_x;M_x,N_x];
            end
        end
    end
    aa =10;
end
