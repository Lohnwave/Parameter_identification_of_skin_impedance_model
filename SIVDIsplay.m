% ����1���������ߡ�2���迹ģ������ʵ�������
% ����ͼ��ͼ����
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
    plot(GbsetData,'k','LineWidth',2) % �߿�Ϊ2
    hold on
    text(10,GbsetData(1)+300,num2str(GbsetData(1)),'color','k');
    hold on
    h=figure(i);
    %axis;
    %axis tight;                     % ��������
    %axis([[0 150000] [0  0.1]]);
    h_axis=get(h,'Children');
    set(h_axis,'LineWidth',1.5);
    
    set(gca,'FontSize',12, 'FontName','Times New Roman');  %��������
    set(gcf,'color','w');   % ������ɫ
%     xlabel('Iteration ','fontsize',16);  %x����  ����12 
    xlabel('���� ','fontsize',16,'FontName','����');  %x����  ����12   
    ylabel(['\fontname{����}���',' / \fontname{Times New Roman}��'],'fontsize',16);
%     ylabel(['Error','\rm/ ��'],'fontsize',16);  
%     title('Convergence curve','fontsize',20,'FontName','Times New Roman');  
    title('��������','fontsize',16,'FontName','����');  
    hl=legend(Fig_algorithm_Name);   %���Ͻǵı�
    %      set(hl,'Box','off');
    % �ֲ��Ŵ� ������ͼ
    %                ���½����      ��ȸ߶�
    axes('Position',[0.48,0.25,     0.28,0.25]); % ������ͼ  
    plot(800:1:1000,GbsetData(800:1000),'k','LineWidth',2);
    hold on
    text(1000,GbsetData(1000),num2str(GbsetData(1000)),'color','k');
    
%     axis([[600 1000] [lb(i) ub(i)]]);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.49 3.37]);
    print('-dtiff','-r600',['SIVConvergence curve����','.tiff']);
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
set(gca,'FontSize',12, 'FontName','Times New Roman');  %��������
set(gcf,'color','w');   % ������ɫ
%     xlabel('f ','fontsize',16);  %x����  ����12 
%     xlabel('\it\fontname{Times New Roman}f \rm(\it\fontname{Times New
%     Roman}x\rm)') % f(x)
%     ������ ͬʱ��ʾб�������
xlabel(['\itf ','\rm/ Hz'],'fontsize',16);
% ylabel(['Modulus of impedance','\rm/ ��'],'fontsize',16);
ylabel(['\fontname{����}�迹ģֵ',' / \fontname{Times New Roman}��'],'fontsize',16);
    %ylabel('log(Average Function Value Error)','fontsize',16);
title(['\fontname{Times New Roman}NSEI','\fontname{����}ģ����������ݱȽ�'],'fontsize',16);  
%     hl=legend('test');   %���Ͻǵı�   
h1=legend('\fontname{����}��������',['\fontname{Times New Roman}NSEI','\fontname{����}ģ��']);
set(h1,'FontSize',12,'FontWeight','normal')
% �ֲ��Ŵ� ������ͼ
%                ���½����      ��ȸ߶�
axes('Position',[0.60,0.25,     0.28,0.25]); % ������ͼ  
plot(f(52:79),Z(52:79),'r*','LineWidth',1);
hold on 
plot(f(52:79),Z0(52:79),'b','LineWidth',2);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.49 3.37]);
print('-dtiff','-r600',['SIVModelVsData����','.tiff']);
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
set(gca,'FontSize',12, 'FontName','Times New Roman');  %��������
set(gcf,'color','w');   % ������ɫ
%     xlabel('f ','fontsize',16);  %x����  ����12 
%     xlabel('\it\fontname{Times New Roman}f \rm(\it\fontname{Times New
%     Roman}x\rm)') % f(x)
%     ������ ͬʱ��ʾб�������
xlabel(['\itf ','\rm/ Hz'],'fontsize',16);
ylabel(['\fontname{����}��ֵ',' / \fontname{Times New Roman}��'],'fontsize',16);
% ylabel(['\fontname{����}����ֵ \fontname{Euclid}this is english!\fontname{����}�������ģ� '])
% Rp = a + b.*exp(c.*(f).^-1);
title(['NSEI: ','\itRp','=','\ita','+','\itbe','^{c/f}'],'fontsize',16,'FontName','Times New Roman');  
%     hl=legend('test');   %���Ͻǵı�   
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.49 3.37]);
print('-dtiff','-r600',['SIVRp����','.tiff']);
