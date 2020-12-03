clear
clc
set(groot,'defaulttextinterpreter','latex');

my_total_data_force=xlsread('data1.csv');
data_1=xlsread('data1.csv');
data_2=xlsread('data2.csv');
data_3=xlsread('data3.csv');
[wide,len]=size(data_1);

my_gamma=data_1(1,2:len);
% my_alphaz=data_1(2:2+6,2:len);
my_beta=data_1(2:2+6,1);
alphaz_1=data_1(2:2+6,2:len);
alphaz_2=data_2(2:2+6,2:len);
alphaz_3=data_3(2:2+6,2:len);
alphax_1=data_1(10:10+6,2:len);
alphax_2=data_2(10:10+6,2:len);
alphax_3=data_3(10:10+6,2:len);
flip_gamma = [fliplr(-pi-my_gamma(:,2:7)) my_gamma fliplr(pi-my_gamma(:,8:13))];
gamma_interval = (flip_gamma(1,end)-flip_gamma(1,1))/25;
gamma_range = flip_gamma(1,1):gamma_interval:flip_gamma(1,end);

alphaz_my=(data_1(2:2+6,2:len)+data_2(2:2+6,2:len)+data_3(2:2+6,2:len))/3;
alphax_my=(data_1(10:10+6,2:len)+data_2(10:10+6,2:len)+data_3(10:10+6,2:len))/3;
% define the full alphaz
flip_alphaz_my = [fliplr(alphaz_my(:,2:7)) alphaz_my fliplr(alphaz_my(:,8:13))];
flip_alphax_my = [-fliplr(alphax_my(:,2:7)) alphax_my -fliplr(alphax_my(:,8:13))];
matrix_size=size(alphaz_my);
%,matrix_size(1)*2,matrix_size(2)*2

total_alphaz=[];
total_alphax=[];
total_gamma=[];
total_beta=[];
[row_n col_n]=size(alphaz_my);
F_gamma= [];
F_beta = [];

fit_gamma = gamma_range;
%-pi:2*pi/25:pi;
fit_num = 0: 25/25 : 25;
% fit_gg = fit (fit_num',fit_gamma','poly1')
fit_gg = fit (fit_gamma',fit_num','poly1')

fit_beta = pi/2 : -pi/6 : -pi/2;
fit_num_b = 0 : 6/6 : 6;
fit_beta = fit(fit_beta',fit_num_b','poly1')
% fit_gg
% plot(fit_gamma,fit_gg(fit_gamma))
% hold on 
% plot(fit_gamma,fit_num,'o')
left_range =  2 : 7;
right_range = 8 : 13;
rep_gamma = [fliplr(-pi-my_gamma(left_range)) my_gamma fliplr(pi-my_gamma(right_range))];
for i = 1 : row_n
%%%%%%% Flip the data to form a periodic dataset
%     rep_gamma = [my_gamma(1:7)-pi/2 my_gamma my_gamma(8:14)+pi/2];
    

%     total_alphaz=[total_alphaz flip(alphaz_1(i,1:7)) alphaz_1(i,:) flip(alphaz_1(i,8:14)) flip(alphaz_2(i,1:7)) alphaz_2(i,:) flip(alphaz_2(i,8:14)) flip(alphaz_3(i,1:7)) alphaz_3(i,:) flip(alphaz_3(i,8:14))];
%     total_alphax=[total_alphax -flip(alphax_1(i,1:7)) alphax_1(i,:) -flip(alphax_1(i,8:14)) -flip(alphax_2(i,1:7)) alphax_2(i,:) -flip(alphax_2(i,8:14)) -flip(alphax_3(i,1:7)) alphax_3(i,:) -flip(alphax_3(i,8:14))];
    total_alphaz=[total_alphaz flip(alphaz_1(i,left_range)) alphaz_1(i,:) flip(alphaz_1(i,right_range)) flip(alphaz_2(i,left_range)) alphaz_2(i,:) flip(alphaz_2(i,right_range)) flip(alphaz_3(i,left_range)) alphaz_3(i,:) flip(alphaz_3(i,right_range))];
    total_alphax=[total_alphax -flip(alphax_1(i,left_range)) alphax_1(i,:) -flip(alphax_1(i,right_range)) -flip(alphax_2(i,left_range)) alphax_2(i,:) -flip(alphax_2(i,right_range)) -flip(alphax_3(i,left_range)) alphax_3(i,:) -flip(alphax_3(i,right_range))];

%     total_alphax=[total_alphax -flip(alphax_1(i,1:7)) alphax_1(i,:) -flip(alphax_1(i,8:14)) -flip(alphax_2(i,1:7)) alphax_2(i,:) -flip(alphax_2(i,8:14)) -flip(alphax_3(i,1:7)) alphax_3(i,:) -flip(alphax_3(i,8:14))];
%     total_alphax=[total_alphax flip(alphax_1(i,1:7)) alphax_1(i,:) flip(alphax_1(i,8:14)) flip(alphax_2(i,1:7)) alphax_2(i,:) flip(alphax_2(i,8:14)) flip(alphax_3(i,1:7)) alphax_3(i,:) flip(alphax_3(i,8:14))];
%     total_gamma=[total_gamma my_gamma(1:7)-pi/2 my_gamma my_gamma(8:14)+pi/2 my_gamma(1:7)-pi/2 my_gamma my_gamma(8:14)+pi/2 my_gamma(1:6)-pi/2 my_gamma my_gamma(8:14)+pi/2];
%%%% gamma and beta for one row   
    total_gamma =[total_gamma repmat(rep_gamma,1,3)];
    total_beta=[total_beta ones(1,length(rep_gamma)*3)*my_beta(i)];
%     F_gamma = [F_gamma; my_gamma(1:7)-pi/2 my_gamma my_gamma(8:14)+pi/2];
    F_gamma = [F_gamma; rep_gamma];
    F_beta = [F_beta; ones(1,length(rep_gamma))*my_beta(i)];
%%%%%%%The original data set
%     total_alphaz=[total_alphaz alphaz_1(i,:) alphaz_2(i,:) alphaz_3(i,:)];
%     total_alphax=[total_alphax alphax_1(i,:) alphax_2(i,:) alphax_3(i,:)];
%     total_gamma=[total_gamma  my_gamma  my_gamma  my_gamma ];
%     total_beta=[total_beta ones(1,col_n*3)*my_beta(i)];
end
%%%%% This part conduct FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change the fitting order here%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
order =1;
[az_co,ax_co,mn_alphaz,mn_alphax,MN_shift_z,MN_shift_x,filtaz,filtax]=JuntaoFFT(order);

% [F_az , F_ax]=JuntaoFFT(1)
%%%% This part gets the Fourier fitting alphaz
M = 7;
N = 26;%28;
% [F_gamma,F_beta]=meshgrid(-pi:0.1:pi,-pi/2:0.1:pi/2); 


F_az = zeros(7,26); %(7,28)
F_ax = zeros(7,26); %(7,28)
coe_z=zeros(length(az_co),5);
coe_x=zeros(length(ax_co),5);
% [Test_gamma,Test_beta] = meshgrid(-pi:0.1:pi,pi/2:-0.1:-pi/2);
[Test_gamma,Test_beta] = meshgrid(rep_gamma,pi/2:-pi/6:-pi/2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Calculate the coefficients for the fitting formula%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:length(az_co)
%     m = mn_alphaz(i,1) -1;
%     n = mn_alphaz(i,2) -1;
    m = MN_shift_z(i,1);
    n = MN_shift_z(i,2);
    
%     coe_z(i,:) = [-pi*2i*m/M*6/pi,pi*2i*n/N*27/(2*pi),az_co(i)*exp(pi*2i*m/M*3+pi*2i*n/N*27/2),m,n];


%     coe_z(i,:) = [-pi*2i*m/M*6/pi,pi*2i*n/N*4.39,az_co(i)*exp(pi*2i*m/M*3+pi*2i*n/N*27/2),m,n];
%     coe_z(i,:) = [-pi*2i*m/M*6/pi,pi*2i*n/N*3.979,az_co(i)*exp(pi*2i*m/M*3+pi*2i*n/N*12.5),m,n];
    coe_z(i,:) = [-pi*2i*m/M*7/pi,pi*2i*n/N*4.249,az_co(i)*exp(pi*2i*m/M*3+pi*2i*n/N*12.5),m,n];
    
    
    %         coe_z(i,:) = [-pi*2i*m/M*6/pi,pi*2i*n/N*4.414,az_co(i)*exp(pi*2i*m/M*3+pi*2i*n/N*12.5),m,n];

%     coe_z(i,:) = [-2i*m,1i*n,az_co(i)*exp(1i*m*pi+pi*1i*n),m,n];
    
end

 for i = 1: length(ax_co)
%     m = mn_alphax(i,1) -1;
%     n = mn_alphax(i,2) -1;
    m = MN_shift_x(i,1);
    n = MN_shift_x(i,2);
%     coe_x(i,:) = [-pi*2i*m/M*6/pi,pi*2i*n/N*27/(2*pi),ax_co(i)*exp(pi*2i*m/M*3+pi*2i*n/N*27/2),m,n];
%     coe_x(i,:) = [-pi*2i*m/M*6/pi,pi*2i*n/N*4.39,ax_co(i)*exp(pi*2i*m/M*3+pi*2i*n/N*27/2),m,n];
%     coe_x(i,:) = [-pi*2i*m/M*6/pi,pi*2i*n/N*3.979,ax_co(i)*exp(pi*2i*m/M*3+pi*2i*n/N*12.5),m,n];
    coe_x(i,:) = [-pi*2i*m/M*7/pi,pi*2i*n/N*4.249,ax_co(i)*exp(pi*2i*m/M*3+pi*2i*n/N*12.5),m,n];
%     coe_x(i,:) = [-pi*2i*m/M*6/pi,pi*2i*n/N*4.414,ax_co(i)*exp(pi*2i*m/M*3+pi*2i*n/N*12.5),m,n];


%     coe_x(i,:) = [-2i*m,1i*n,ax_co(i)*exp(1i*m*pi+pi*1i*n),m,n];
 end
%%%%%%%%%%Get the coefficients %%%%%%%%%%%%%%%%%%
AB_list = coe_z(:,3);
A_list = real(coe_z(:,3));
% A_list (abs(A_list) < 0.005) = 0;
B_list = -imag(coe_z(:,3));
% B_list (abs(B_list)< 0.005) = 0;
AB_mn = imag(coe_z(:,1:2));

C_list = real(coe_x(:,3));
% C_list (abs(C_list)<0.008) = 0;
D_list = -imag(coe_x(:,3));
% D_list (abs(D_list)<0.005) = 0;
CD_mn = imag(coe_x(:,1:2));

alphazz = 0;
alphaxx = 0;
fit_z = 0;
fit_x = 0;
Fit_x = 0;
for i = 1:length(az_co)
    
%     a_b = round(AB_mn(i,1))*Test_beta+round(AB_mn(i,2))*Test_gamma;
    a_b = AB_mn(i,1)*Test_beta+AB_mn(i,2)*Test_gamma;
    a_b_error = AB_mn(i,1)*total_beta+AB_mn(i,2)*total_gamma;
    
    alphazz = alphazz + A_list(i)*cos(a_b) + B_list(i) * sin(a_b); 
    fit_z = fit_z + A_list(i)*cos(a_b_error) + B_list(i) * sin(a_b_error); 
%     alphazz = alphazz + real(AB_list(i)*exp(a_b));
end


for i = 1:length(ax_co)
    
%     a_b = round(CD_mn(i,1))*Test_beta+round(CD_mn(i,2))*Test_gamma;
    c_d = CD_mn(i,1)*Test_beta+CD_mn(i,2)*Test_gamma; 
    c_d_error = CD_mn(i,1)*total_beta+CD_mn(i,2)*total_gamma;

    alphaxx = alphaxx + C_list(i)*cos(c_d) + D_list(i) * sin(c_d); 
    fit_x = fit_x + C_list(i)*cos(c_d_error) + D_list(i) * sin(c_d_error);
    Fit_x = Fit_x + C_list(i)*cos(c_d) + D_list(i) * sin(c_d);
%     alphazz = alphazz + real(AB_list(i)*exp(a_b));
end
error_z = sum(abs((fit_z-total_alphaz)))/length(fit_z)
error_x = sum(abs((fit_x-total_alphax)))/length(fit_x)
% error_z = abs(sum((fit_z-total_alphaz)))/length(fit_z)
% error_x = abs(sum((fit_x-total_alphax)))/length(fit_x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Plot the surface with integer index%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1: row_n
    for l = 1:col_n*2-2
        for i = 1:length(az_co)
            Amn = real(az_co(i));
            Bmn = imag(az_co(i));
            m = mn_alphaz(i,1) -1;
            n = mn_alphaz(i,2) -1;


            beta_current = F_beta(k,l);
            gamma_current = F_gamma(k,l);
            F_az(k,l) = real(F_az(k,l) + az_co(i)*exp(pi*2i*m/M*(k-1)+pi*2i*n/N*(l-1)));
%             F_az(k,l) = real(F_az(k,l) + az_co(i)*exp(pi*2i*m/M*(-6/pi*beta_current+3)+pi*2i*n/N*(3.979*gamma_current+12.5)));

%             F_az(k,l) = real(F_az(k,l) + az_co(i)*exp(pi*2i*m/M*(-6/pi*beta_current+3)+pi*2i*n/N*(27/(2*pi)*gamma_current+27/2)));
%             F_az(k,l) = real(F_az(k,l) + az_co(i)*exp(pi*2i*m/M*(-6/pi*beta_current+3)+pi*2i*n/N*(4.39*gamma_current+27/2)));

        end
        
        for i = 1: length(ax_co)
            Cmn = real(ax_co(i));
            Dmn = imag(ax_co(i));
            m = mn_alphax(i,1) -1;
            n = mn_alphax(i,2) -1;

            F_ax(k,l) = real(F_ax(k,l) + ax_co(i)*exp(pi*2i*m/M*(k-1)+pi*2i*n/N*(l-1)));    
%             F_ax(k,l) = real(F_ax(k,l) + ax_co(i)*exp(pi*2i*m/M*(-6/pi*beta_current+3)+pi*2i*n/N*(3.979*gamma_current+12.5)));       
%             F_ax(k,l) = real(F_ax(k,l) + ax_co(i)*exp(pi*2i*m/M*(-6/pi*beta_current+3)+pi*2i*n/N*(27/(2*pi)*gamma_current+27/2)));
            


%             F_ax(k,l) = real(F_ax(k,l) + ax_co(i)*exp(2i*m*beta_current+1i*n*gamma_current));
%              F_ax(k,l) = real(F_ax(k,l) + ax_co(i)*exp(pi*2i*m/M*(-6/pi*beta_current+3)+pi*2i*n/N*(4.39*gamma_current+27/2)));
        end
%          fprintf('al=%f\n',27/(2*pi)*gamma_current+27/2)
    end
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  for i = 1:length(az_co)
%     Amn = real(az_co(i));
%     Bmn = imag(az_co(i));
%     m = mn_alphaz(i,1) -1;
%     n = mn_alphaz(i,2) -1;
%     F_az = real(F_az + az_co(i)*exp(pi*2i*m/32*(6/pi*F_beta+3)+pi*2i*n/63*(27/(2*pi)*F_gamma+27/2)));
% 
% end
% 
% for i = 1: length(ax_co)
%     Cmn = real(ax_co(i));
%     Dmn = imag(ax_co(i));
%     m = mn_alphax(i,1) -1;
%     n = mn_alphax(i,2) -1;
% %     F_ax = real(F_ax + ax_co(i)*exp(pi*2i*m/M*(k-1)+pi*2i*n/N*(l-1)));    
% end
% F_az = az_co(order+1,order+1);
% F_ax = ax_co(order+1,order+1);
% for j = 1: length(az_co)
%     for i = 1: length(az_co)
%         Amn = real(az_co(i,j));
%         Bmn = imag(az_co(i,j));
%         Cmn = real(ax_co(i,j));
%         Dmn = imag(ax_co(i,j));
%         m = order+1-i;
%         n = j-order-1;
%         
%         fprintf('m=%f, n=%f,  Amn=%f\n',m,n,Amn)
%         F_az = F_az + Amn* cos(2*m*F_beta+n*F_gamma)+ Bmn * sin(2*m*F_beta+n*F_gamma);
%         F_ax = F_ax + Cmn* cos(2*m*F_beta+n*F_gamma)+ Dmn * sin(2*m*F_beta+n*F_gamma);
%         
%     end  
% end

% surface(F_beta,F_gamma,F_az)
% [fit_alpha_z good_z]=fit([total_beta' total_gamma'],total_alphaz','poly23');
% [fit_alpha_x good_x]=fit([total_beta' total_gamma'],total_alphax','poly23');
% 
% 
% alphax_row=[];
% alphaz_row=[];
% poly_fitted_ax=ones(row_n,col_n);
% poly_fitted_az=ones(row_n,col_n);
% Fourier_fitted_ax=ones(row_n,col_n);
% Fourier_fitted_az=ones(row_n,col_n);
% for i = 1 : row_n
% 
%     for j= 1: col_n
%         poly_fitted_ax(i,j)=fit_alpha_x(my_beta(i),my_gamma(j));
%         poly_fitted_az(i,j)=fit_alpha_z(my_beta(i),my_gamma(j));
%         Fourier_fitted_ax(i,j)=F_fitted_alphax(my_beta(i),my_gamma(j));
%         Fourier_fitted_az(i,j)=F_fitted_alphaz(my_beta(i),my_gamma(j));
% %         fprintf('beta= %f, gamma= %f \n',my_beta(i),my_gamma(j))
%     end
%     
% end

% poly_error_x=poly_fitted_ax-alphax_my;
% poly_error_z=poly_fitted_az-alphaz_my;
% Fourier_error_x=Fourier_fitted_ax-alphax_my;
% Fourier_error_z=Fourier_fitted_ax-alphax_my;


%%%%%% This part plot the color map of alpha x and alpha z for each beta
%%%%%% and gamma pair
% figure
% subplot(1,2,1)
% xvalues=categorical(round(my_gamma*180/pi));
% yvalues=categorical(round(my_beta.'*180/pi));
% 
% hx=heatmap(xvalues,yvalues,Fourier_fitted_az,'Colormap',jet);
% 
% hx.Title = {'\alpha_z'};
% 
% hx.XLabel = '\gamma';
% hx.YLabel = '\beta';
% hx.FontSize=20;
% caxis([-0.25,0.25]);
% 
% subplot(1,2,2)
% hz=heatmap(xvalues,yvalues,Fourier_fitted_ax,'Colormap',jet);
% hz.Title = '\alpha_x';
% 
% hz.XLabel = '\gamma';
% hz.YLabel = '\beta';
% caxis([-0.11,0.11]);
% 
% hz.FontSize=20;


%%%%%%%%%%%% This part plots the colormap errors %%%%%%%%%%%%%%%
% figure
% 
% subplot(1,2,1)
% xvalues=categorical(my_gamma);
% yvalues=categorical(my_beta.');
% 
% hx=heatmap(xvalues,yvalues,poly_error_z,'Colormap',jet);
% 
% hx.Title = {'Error \alpha_z'};
% 
% hx.XLabel = '\gamma';
% hx.YLabel = '\beta';
% hx.FontSize=20;
% caxis([-0.1,0.1]);
% % caxis([-0.25 0.25]);
% subplot(1,2,2)
% hz=heatmap(xvalues,yvalues,poly_error_x,'Colormap',jet);
% hz.Title = 'Error \alpha_x';
% 
% hz.XLabel = '\gamma';
% hz.YLabel = '\beta';
% caxis([-0.1,0.1]);
% 
% hz.FontSize=20;
% 
% 
% colormap([[(0:127)/127; (0:127)/127; ones(1,128)],[ones(1,128); (127:-1:0)/127; (127:-1:0)/127]]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%% This part plots raw data and fitted surface comparisions  %%%
figure 

% plot3(total_beta,total_gamma,total_alphaz,'o','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r')
% hold on
% plot3(total_beta,total_gamma,fit_z,'+','MarkerSize',5)
plot3(total_beta,total_gamma,fit_z-total_alphaz,'o','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r')
% plot(fit_alpha_z,[total_beta' total_gamma'],total_alphaz')
% title('Fitted $\alpha_z$','FontSize',30)
% plot3(F_beta,F_gamma,F_az,'+','MarkerSize',5)
% plot3(Test_beta,Test_gamma,alphazz,'+','MarkerSize',5)
% surf(Test_beta,Test_gamma,alphazz,'EdgeAlpha',0)
% surf(F_beta,F_gamma,F_az,'EdgeAlpha',0)
surf(F_beta,F_gamma,F_az-flip_alphaz_my,'EdgeAlpha',0)
alpha(0.4)
xlabel('$\beta$','FontSize',20)
ylabel('$\gamma$','FontSize',20)
zlabel('$\alpha_z$','FontSize',20)
% title(['$\alpha_z $ Fitting order:',num2str(order)])
title(['$\alpha_z $ errors order:',num2str(order)])
set(gca,'FontSize',20)
colormap([[(0:127)/127; (0:127)/127; ones(1,128)],[ones(1,128); (127:-1:0)/127; (127:-1:0)/127]]')
figure 
% plot3(total_beta,total_gamma,total_alphax,'o','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r')
% hold on
% plot3(total_beta,total_gamma,fit_x,'+','MarkerSize',5)
plot3(total_beta,total_gamma,fit_x-total_alphax,'o','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r')
% plot3(F_beta,F_gamma,F_ax,'+','MarkerSize',5)
% surf(Test_beta,Test_gamma,alphaxx,'EdgeAlpha',0)
% surf(F_beta,F_gamma,F_ax,'EdgeAlpha',0)
surf(F_beta,F_gamma,Fit_x-flip_alphax_my,'EdgeAlpha',0)
alpha(0.4)
% plot(fit_alpha_x,[total_beta' total_gamma'],total_alphax')
% title('Fitted $\alpha_x$','FontSize',30)
xlabel('$\beta$','FontSize',20)
ylabel('$\gamma$','FontSize',20)
zlabel('$\alpha_x$','FontSize',20)
% title(['$\alpha_x $ Fitting order:',num2str(order)])
title(['$\alpha_x $ errors order:',num2str(order)])
set(gca,'FontSize',20)
colormap([[(0:127)/127; (0:127)/127; ones(1,128)],[ones(1,128); (127:-1:0)/127; (127:-1:0)/127]]')
%%%%%%%%%%%%%% This part plot the fitting data errors%%%%%%%%%%%%%%%%%%%
figure

subplot(1,2,1)

imagesc([-pi,pi],[-pi/2,pi/2],(Fit_x-flip_alphax_my)/mean(abs(total_alphax)));
axis equal
axis([-pi,pi,-pi/2,pi/2]);
colorbar
% caxis([-0.1,0.1]);
xlabel('$\gamma$')
ylabel('$\beta$');
title(['$\alpha_x$ error Order:',num2str(order)]);
set(gca,'FontSize',20)
subplot(1,2,2)

imagesc([-pi,pi],[-pi/2,pi/2],(F_az-flip_alphaz_my)/mean(abs(total_alphaz)));
axis equal
axis([-pi,pi,-pi/2,pi/2]);
colorbar
% caxis([-0.1,0.1]);
xlabel('$\gamma$')
ylabel('$\beta$');
title(['$\alpha_z$ error Order:',num2str(order)]);
set(gca,'FontSize',20)
colormap([[(0:127)/127; (0:127)/127; ones(1,128)],[ones(1,128); (127:-1:0)/127; (127:-1:0)/127]]')



% colormap([[(0:127)/127; (0:127)/127; ones(1,128)],[ones(1,128); (127:-1:0)/127; (127:-1:0)/127]]')
