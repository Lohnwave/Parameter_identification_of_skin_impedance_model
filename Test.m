
clear
clc;

% mex  cec13_func.cpp -DWINDOWS 
D=7;
M=10;
% Xmin=-100;
% Xmax=100;
% Xmin = [0, 35.39*10^-6, 35.39*10^-6, -10000, -10000, -10000, 0];
Xmin = [0, 0, 0, -10000, -10000, -10000, 0];
% Xmax = [4500, 0.1592, 0.1592, 4500, 4500, 4500, 4500];
Xmax = [4500, 0.1592, 0.1592, 10000, 10000, 10000, 4500];
% Xmax = [10000, 0.1592, 0.1592, 10000, 10000, 10000, 10000];
pop_size=50;
iter_max=1000;
max_fes=50000;
warning off;
runs=2;

algorithm_Name={'EBS_SCA_ZK'};
% fhd=str2func('cec13_func');
Data1=zeros(runs,size(algorithm_Name,2));   % Store Error
Data2=zeros(runs,size(algorithm_Name,2));   % Store the Fitness value of the Global optimum solution
Data3=zeros(runs,iter_max,size(algorithm_Name,2));      % Store the fitness value of the best individual in the search process
Data4=zeros(runs,D,size(algorithm_Name,2)); % Store the Global optimum solution
t=zeros(runs,size(algorithm_Name,2));
f_mean=zeros(size(algorithm_Name,2));
t_mean=zeros(size(algorithm_Name,2));
E_mean=zeros(size(algorithm_Name,2));
for k=1:size(algorithm_Name,2)
    fprintf('Algorithm =\t %d\n',k);
        for j=1:runs
            fprintf('run =\t %d\n',j);
            BSOFUNC=algorithm_Name{k};
            tic;
            [gbestval,allgbestval,gbest,fitcount]=feval(BSOFUNC,iter_max,max_fes,pop_size,D,Xmin,Xmax);
            t(j,k)=toc;
%             Data1(i,j,k)=Error;
            Data2(j,k)=gbestval;
            Data3(j,:,k)=allgbestval;
            Data4(j,:,k)=gbest;
        end
        f_mean(k)=mean(Data2(:,k));
        t_mean(k)=mean(t(:,k));
        E_mean(k)=mean(Data1(:,k));
    file_name= [algorithm_Name{k},'_M10_7D_',num2str(runs),'runs.mat'];
    save (file_name);
end
load('skin_impedanceData.mat')
%%%
