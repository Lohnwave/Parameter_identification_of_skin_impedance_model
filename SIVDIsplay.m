% 绘制1、收敛曲线。2、阻抗模型与真实数据误差
% 包含图中图绘制
clear;
close all;
%% Convergence curve
load('SIV_ZK_M10_7D_30runs.mat')
Fig_algorithm_Name=char('EBS-SCA');
%GbsetData=zeros(28,iter_max,size(algorithm_Name,1));
% GbsetData = mean(Data3,1);
[minError,GbsetDataindex] = min(Data3(:,1000));
GbsetData = Data3(GbsetDataindex,:);
for i=1:1
    figure(i) 
    plot(GbsetData,'k','LineWidth',2) % 线宽为2
    hold on
    text(10,GbsetData(1)+300,num2str(GbsetData(1)),'color','k');
    hold on
    h=figure(i);
    %axis;
    %axis tight;                     % 紧坐标轴
    %axis([[0 150000] [0  0.1]]);
    h_axis=get(h,'Children');
    set(h_axis,'LineWidth',1.5);
    
    set(gca,'FontSize',12, 'FontName','Times New Roman');  %设置字体
    set(gcf,'color','w');   % 背景白色
%     xlabel('Iteration ','fontsize',16);  %x坐标  字体12 
    xlabel('迭代 ','fontsize',16,'FontName','宋体');  %x坐标  字体12   
    ylabel(['\fontname{宋体}误差',' / \fontname{Times New Roman}Ω'],'fontsize',16);
%     ylabel(['Error','\rm/ Ω'],'fontsize',16);  
%     title('Convergence curve','fontsize',20,'FontName','Times New Roman');  
    title('收敛曲线','fontsize',16,'FontName','宋体');  
    hl=legend(Fig_algorithm_Name);   %右上角的标
    %      set(hl,'Box','off');
    % 局部放大 绘制子图
    %                左下角左边      宽度高度
    axes('Position',[0.48,0.25,     0.28,0.25]); % 生成子图  
    plot(800:1:1000,GbsetData(800:1000),'k','LineWidth',2);
    hold on
    text(1000,GbsetData(1000),num2str(GbsetData(1000)),'color','k');
    
%     axis([[600 1000] [lb(i) ub(i)]]);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.49 3.37]);
    print('-dtiff','-r600',['SIVConvergence curve中文','.tiff']);
end
%% NSEI model Vs real data
load('skin_impedanceData.mat')
figure(2)
plot(f,Z,'r*','LineWidth',1);
hold on
Z0 = PLOTfobj(Data4(GbsetDataindex,:));
plot(f,Z0,'b','LineWidth',2);
% grid on
hold on 

h=figure(2);
h_axis=get(h,'Children');
set(h_axis,'LineWidth',1.5);
set(gca,'FontSize',12, 'FontName','Times New Roman');  %设置字体
set(gcf,'color','w');   % 背景白色
%     xlabel('f ','fontsize',16);  %x坐标  字体12 
%     xlabel('\it\fontname{Times New Roman}f \rm(\it\fontname{Times New
%     Roman}x\rm)') % f(x)
%     横坐标 同时显示斜体和正体
xlabel(['\itf ','\rm/ Hz'],'fontsize',16);
% ylabel(['Modulus of impedance','\rm/ Ω'],'fontsize',16);
ylabel(['\fontname{宋体}阻抗模值',' / \fontname{Times New Roman}Ω'],'fontsize',16);
    %ylabel('log(Average Function Value Error)','fontsize',16);
title(['\fontname{Times New Roman}NSEI','\fontname{宋体}模型与测量数据比较'],'fontsize',16);  
%     hl=legend('test');   %右上角的标   
h1=legend('\fontname{宋体}测量数据',['\fontname{Times New Roman}NSEI','\fontname{宋体}模型']);
set(h1,'FontSize',12,'FontWeight','normal')
% 局部放大 绘制子图
%                左下角左边      宽度高度
axes('Position',[0.60,0.25,     0.28,0.25]); % 生成子图  
plot(f(52:79),Z(52:79),'r*','LineWidth',1);
hold on 
plot(f(52:79),Z0(52:79),'b','LineWidth',2);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.49 3.37]);
print('-dtiff','-r600',['SIVModelVsData中文','.tiff']);
RMSE = sqrt(mean((Z0-Z).^2));
corr = corrcoef(Z0,Z);
%% Rp curve
% Rp = a + b.*exp(c.*(f).^-1);
a = Data4(GbsetDataindex,4);
b = Data4(GbsetDataindex,5);
c = Data4(GbsetDataindex,6);
Rp = a + b.*exp(c.*(f).^-1);
figure(3)
plot(f,Rp,'k','LineWidth',2);

h=figure(3);
h_axis=get(h,'Children');
set(h_axis,'LineWidth',1.5);
set(gca,'FontSize',12, 'FontName','Times New Roman');  %设置字体
set(gcf,'color','w');   % 背景白色
%     xlabel('f ','fontsize',16);  %x坐标  字体12 
%     xlabel('\it\fontname{Times New Roman}f \rm(\it\fontname{Times New
%     Roman}x\rm)') % f(x)
%     横坐标 同时显示斜体和正体
xlabel(['\itf ','\rm/ Hz'],'fontsize',16);
ylabel(['\fontname{宋体}阻值',' / \fontname{Times New Roman}Ω'],'fontsize',16);
% ylabel(['\fontname{宋体}理论值 \fontname{Euclid}this is english!\fontname{宋体}又是中文！ '])
% Rp = a + b.*exp(c.*(f).^-1);
title(['NSEI: ','\itRp','=','\ita','+','\itbe','^{c/f}'],'fontsize',16,'FontName','Times New Roman');  
%     hl=legend('test');   %右上角的标   
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.49 3.37]);
print('-dtiff','-r600',['SIVRp中文','.tiff']);
