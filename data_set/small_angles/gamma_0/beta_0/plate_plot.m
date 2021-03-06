clc
clear
%%%%%%%% create two arrays to store time move down and time move
%%%%%%%% horizontally
tdown=zeros(1,10);
thorizontal=zeros(1,10);
time_settle=0.5;
intrude_period=0.2;

for k =1:11
    tdown(k)=time_settle+1.2*(k-1);
    thorizontal(k)=time_settle+intrude_period+1.2*(k-1);
end

force_address='output_plate_forces.csv';
position_address='output_plate_positions.csv';
total_data_force=csvread(force_address);
total_data_position=csvread(position_address);
stress_list=[];
% slope_list=[];
for i =1:10
    subplot(4,3,i)
    start_index=find_index(1,thorizontal(i),total_data_force(:,1))+250;
    term_index=find_index(1,tdown(i+1),total_data_force(:,1))-10;
    range=start_index:term_index;
    time=total_data_force(range,1);
    f_x=total_data_force(range,2);
    f_y=total_data_force(range,3);
    f_z=total_data_force(range,4);
    gamma=total_data_force(1,1);
    beta=total_data_force(1,2);
    pos_z=total_data_position(:,4);
    pos_x=total_data_position(:,2);
    gm_top=total_data_position(start_index,5);
    dep=gm_top-pos_z(start_index);
    Area=3.81*2.54;
    stress_current=[dep;sum(f_x)/(term_index-start_index+1)/Area;sum(f_z)/(term_index-start_index+1)/Area];
    if i>3 & i<10
        stress_list=[stress_list,stress_current];
    end
    
    plot(time,f_x/Area);
    hold on
    plot(time,f_y/Area);
    hold on
    plot(time,f_z/Area);
    title("\gamma: "+ gamma + " \beta: "+beta+"  depth:"+dep+ " cm");
    xlabel('time (s)','FontSize',15);
    ylabel('Stress (N/cm^2)','FontSize',15);
    legend('Stress x','Strees y','Stress z','FontSize',10);
    set(gca,'FontSize',15);

end
s_x=polyfit(stress_list(1,:).',stress_list(2,:).',1)
s_z=polyfit(stress_list(1,:).',stress_list(3,:).',1)
figure
plot(total_data_force(:,1),pos_x);
hold on
plot(total_data_force(:,1),pos_z);
hold on
xlabel('time (s)','FontSize',15);
ylabel('m','FontSize',15);
legend('pos x','pos z','FontSize',20);
% start_index=find_index(1,2.5,time);
% for i= 1:6
%     range_start=2+3002*(i-1);
%     range_terminal=3002+3002*(i-1);
%     range=range_start:range_terminal;
%     raw_data_start=1;
%     Area=3.81*2.54;
%     fitted_start=1;
%     fitted_terminal=1;
%     gamma=total_data_force(1+3002*(i-1),1);
%     beta=total_data_force(1+3002*(i-1),2);
%     %create some arrays to store the forces and positions
%     time=total_data_force(range,1);
%     f_x=total_data_force(range,2);
%     f_y=total_data_force(range,3);
%     f_z=total_data_force(range,4);
%     pos_z=total_data_position(range,4);
%     %granular materials top
%     %to find the start index for raw data
%     raw_data_start=find_index(raw_data_start,0.01,f_z);
%     gm_top=total_data_position(raw_data_start+3002*(i-1),4);
%     raw_data_terminal=length(time);
%     raw_data_range=raw_data_start:raw_data_terminal;
%     %to find the terminal index for fitting data
%     fitted_start=find_index(fitted_start,2,gm_top-pos_z);
%     %to find the terminal index for fitting data
%     fitted_terminal=find_index(fitted_terminal,6,gm_top-pos_z);
%     fitted_data_range=fitted_start:fitted_terminal;
% 
%     dep=gm_top-pos_z(fitted_data_range);
%     s_x=polyfit(dep,f_x(fitted_data_range)/Area,1);
%     s_z=polyfit(dep,f_z(fitted_data_range)/Area,1);
%     sx=polyval(s_x,dep);
%     sz=polyval(s_z,dep);
%     %plot figures
%     subplot(3,2,i)
%     plot(gm_top-pos_z(raw_data_range),f_x(raw_data_range)/Area);
%     hold on
%     plot(gm_top-pos_z(raw_data_range),f_z(raw_data_range)/Area);
%     hold on
%     plot(dep,sx,'LineWidth',4);
%     hold on
%     plot(dep,sz,'LineWidth',4);
%     title("$\gamma$: "+ gamma + " $\beta$: "+beta);
%     xlabel('depth (cm)','FontSize',15);
%     ylabel('Stress (N/cm^2)','FontSize',15);
%     legend('Stress x','Strees z','Fitted-stress x','Fitted-stress z','FontSize',10);
%     set(gca,'FontSize',15)
%     slope=[gamma;beta;s_x(1);s_z(1)];
%     slope_list=[slope,slope_list];
% end
% hold off;
% figure
% subplot(1,2,1)
% xvalues=categorical(slope_list(1,:));
% yvalues=categorical(slope_list(2));
% hx=heatmap(xvalues,yvalues,slope_list(3,:),'Colormap',jet);
% hx.Title = '\alpha_x (N/cm^3)';
% hx.XLabel = '\gamma';
% hx.YLabel = '\beta';
% hx.FontSize=20;
% subplot(1,2,2)
% hz=heatmap(xvalues,yvalues,slope_list(4,:),'Colormap',jet);
% hz.Title = '\alpha_z (N/cm^3)';
% hz.XLabel = '\gamma';
% hz.YLabel = '\beta';
% hz.FontSize=20;