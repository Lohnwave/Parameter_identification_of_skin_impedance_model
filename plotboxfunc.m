clear;
close all;
load('SIV_ZK_M10_7D_30runs.mat')
%Note: ***.mat is the result of storing the output (error) of the algorithms 
%      tested on the test set with thirty times.
% data(1,:,:)=EBSMSCA05031;
% data(2,:,:)=GBSO01M10;
% data(3,:,:)=CLPSO;
% data(4,:,:)=ABC;
% data(5,:,:)=ISCA;
% data(6,:,:)=SCA_PSO02;
% B={'EBS-SCA';'GBSO';'CLPSO';'ABC';'ISCA';'SCAPSO'};
Data4 = Data4';
data1(1,:)=Data4(1,:);
data1(2,:)=Data4(7,:);

data2(1,:)=Data4(2,:);%a
data2(2,:)=Data4(3,:);%b

data3(1,:)=Data4(4,:);%a
data3(2,:)=Data4(5,:);%b
data3(3,:)=Data4(6,:);%c



% B={'S0';'MA1';'MA2';'MA3';'MA4';'MA5';'MA6';'MA7';'MA8';'MA9'};
 B1={'Rd','r'};
 B2={'Cd','Cp'};
 B3={'a','b','c'};
for i=1:1
    figure(i)
    A(:,:)=data1(:,:);
    boxplot(A',B1);
    c=boxplot(A',B1);
    set(c,'LineWidth',1.5);
    h=figure(i);
    axis;
    
    h_axis=get(h,'Children');
    set(h_axis,'LineWidth',1.5);
    set(gca,'FontSize',11, 'FontName','Times New Roman');  %设置字体
    set(gcf,'color','w');   % 背景白色
    ylabel('Error','fontsize',16);
    title(['f',num2str(i)],'fontsize',20,'FontName','Times New Roman'); 
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.49 3.37]);
%     print('-dtiff','-r600',[num2str(i),'.tiff']);
end

for i=2:2
    figure(i)
    A2(:,:)=data2(:,:);
    boxplot(A2',B2);
    c=boxplot(A2',B2);
    set(c,'LineWidth',1.5);
    h=figure(i);
    axis;
    
    h_axis=get(h,'Children');
    set(h_axis,'LineWidth',1.5);
    set(gca,'FontSize',11, 'FontName','Times New Roman');  %设置字体
    set(gcf,'color','w');   % 背景白色
    ylabel('Error','fontsize',16);
    title(['f',num2str(i)],'fontsize',20,'FontName','Times New Roman'); 
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.49 3.37]);
%     print('-dtiff','-r600',[num2str(i),'.tiff']);
end

for i=3:3
    figure(i)
    A3(:,:)=data3(:,:);
    boxplot(A3',B3);
    c=boxplot(A3',B3);
    set(c,'LineWidth',1.5);
    h=figure(i);
    axis;
    
    h_axis=get(h,'Children');
    set(h_axis,'LineWidth',1.5);
    set(gca,'FontSize',11, 'FontName','Times New Roman');  %设置字体
    set(gcf,'color','w');   % 背景白色
    ylabel('Error','fontsize',16);
    title(['f',num2str(i)],'fontsize',20,'FontName','Times New Roman'); 
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.49 3.37]);
%     print('-dtiff','-r600',[num2str(i),'.tiff']);
end