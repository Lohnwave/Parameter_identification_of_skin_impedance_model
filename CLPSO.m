%function [gbestval,allgbestval,gbest,fitcount]= CLPSO(fhd,Max_Gen,Max_FES,Particle_Number,Dimension,VRmin,VRmax,varargin)
function [gbestval,allgbestval,gbest,fitcount]= CLPSO(Max_Gen,Max_FES,Particle_Number,Dimension,VRmin,VRmax)
%[gbest,gbestval,fitcount]= CLPSO_new_func('f8',3500,200000,30,30,-5.12,5.12)
%Comprehensive Learning Particle Swarm Optimizer for Global Optimization of
%Multimodal Functions 论文实现代码
rand('state',sum(100*clock));%
me=Max_Gen;%最大迭代次数
ps=Particle_Number;%单次迭代粒子群中粒数
D=Dimension;%单个粒子的空间维度
cc=[1 1];   %acceleration constants
t=0:1/(ps-1):1;t=5.*t;
Pc=0.0+(0.5-0.0).*(exp(t)-exp(t(1)))./(exp(t(ps))-exp(t(1)));
% Pc=0.5.*ones(1,ps);
m=0.*ones(ps,1);
iwt=0.9-(1:me)*(0.7/me);%惯性因子，参见CLPSO论文参考文献[8]
% iwt=0.729-(1:me)*(0.0/me);
cc=[1.49445 1.49445]; %

if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);%粒子空间最小值
    VRmax=repmat(VRmax,1,D);%粒子空间最大值
end

mv=0.2*(VRmax-VRmin);%速度范围限制确定(参数取0-1之间，这里取了0.2)
VRmin=repmat(VRmin,ps,1);%产生1列，ps个VRmin值
VRmax=repmat(VRmax,ps,1);%repmat(VRmax,ps,n)表示产生n列，每列有ps个VRmin
Vmin=repmat(-mv,ps,1);%速度最小值
Vmax=-Vmin;%速度最大值
pos=VRmin+(VRmax-VRmin).*rand(ps,D);%位置初始化

     for i=1:ps;
     %e(i)=feval(fhd,pos(i,:),varargin{:}); % fitness value initialization
        e(i)=fobj(pos(i,:));
     end
%        e=feval(fhd,pos',varargin{:});%D*pop_size matrix.

fitcount=ps; %  evaluation of fitness numbers
vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles
pbest=pos; %pbest initialization
pbestval=e; %initialize the pbest's fitness value
[gbestval,gbestid]=min(pbestval);%gbestval全局最优值，gbestid全局最优值所对应的序号
allgbestval=zeros(1,Max_Gen);
gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
gbestrep=repmat(gbest,ps,1); % sotre gbest for each particle

stay_num=zeros(ps,1);% stagnated iteration

    ai=zeros(ps,D);
    %******************************************
    %1  1  1... 1, D "1"
    %2  2  2... 2, D "2"
    %3  3  3... 3, D "3"
    %|  |  |... |
    %ps ps ps...ps,D "ps"
    f_pbest=1:ps;
    f_pbest=repmat(f_pbest',1,D);% store the index of the excemple for ps particles, and index from the row of the f_pbest represent the corrsponding dimension is learned from that of selected particle. 
    %*******************************************
        for k=1:ps       
            ar=randperm(D);% generate randomly D numbers ranging from 1 to D
            ai(k,ar(1:m(k)))=1;% each number is zeros in ai.
        %************************************
        %using the martix not loop to improve the computing efficiency
        %determining which dimension of particles'pbest selcted randomly the updated particle's each dimension should follow. 
        % updated particle i=(pbesti1,pbesti2,...,pbestiD)
        %fi1=(2,3,....5);fi2(D,1,...,4)
        %pbesti1 is updated using the first dimension of the winner between particle 2 and particle D
        %pbesti2 is updated using the second dimension of the winner between particle 3 and particle 1 
            fi1=ceil(ps*rand(1,D));% select randomly D particles for the D dimension of the updated particle
            fi2=ceil(ps*rand(1,D));% reselect randomly D particles for the D dimension of the updated particle
            fi=(pbestval(fi1)<pbestval(fi2)).*fi1+(pbestval(fi1)>=pbestval(fi2)).*fi2; % select the corrsponding dimension of the particle owning the higher fitness as the excemple to learn form for that dimension.
        %************************************
            bi=ceil(rand(1,D)-1+Pc(k));% if rand()<pc is ture, the corresponding dimension of bi =1, or else is 0
            if bi==zeros(1,D),rc=randperm(D);bi(rc(1))=1;end % if bi is zeros, this means that all exemplars of a particle are its own pbest, we will randomly choose one dimension to learn form antoher particle's pbest's corrsponding dimension
            f_pbest(k,:)=bi.*fi+(1-bi).*f_pbest(k,:);% creat the index of the ps excemple for the particle updated. each row represent the updated particle and each column is the index of the selected particle to replace the updated the corresponding dimension of the particle .
         %************************************
        end

    i=1;


while i<=me&fitcount<=Max_FES %corresponding 129 line
     i=i+1;
            for k=1:ps % corresponding 108 line

                    if stay_num(k)>=5
                %     if round(i/10)==i/10%|stay_num(k)>=5
                        stay_num(k)=0;
                        ai(k,:)=zeros(1,D);
                        f_pbest(k,:)=k.*ones(1,D); 
                        ar=randperm(D);
                        ai(k,ar(1:m(k)))=1;
                        fi1=ceil(ps*rand(1,D));
                        fi2=ceil(ps*rand(1,D));
                        fi=(pbestval(fi1)<pbestval(fi2)).*fi1+(pbestval(fi1)>=pbestval(fi2)).*fi2;
                        bi=ceil(rand(1,D)-1+Pc(k));
                        if bi==zeros(1,D)
                            rc=randperm(D);
                            bi(rc(1))=1;
                        end
                        f_pbest(k,:)=bi.*fi+(1-bi).*f_pbest(k,:);%
                    end

                    for dimcnt=1:D
                        %*******************
                         %******************
                        %1  1  1... 1  D "1" represent the excemple for the first updated particle
                        %2  2  2... 2, D "2" represent the excemple for the second updated particle
                        %3  3  3... 3, D "3"
                        %|  |  |... |
                        %ps ps ps...ps,D "ps"
                        if f_pbest(k,dimcnt)==0
                            f_pbest(k,dimcnt)=ceil(ps*rand());
                        end
                        pbest_f(k,dimcnt)=pbest(f_pbest(k,dimcnt),dimcnt);%f_pbest(k,dimcnt) means which pbest is used as the excemple to update the dimcnt dimension of the particle. 
                    end
                    
                   % aa(k,:)=cc(1).*(1-ai(k,:)).*rand(1,D).*(pbest_f(k,:)-pos(k,:))+cc(2).*ai(k,:).*rand(1,D).*(gbestrep(k,:)-pos(k,:));%~~
                   aa(k,:)=cc(1).*rand(1,D).*(pbest_f(k,:)-pos(k,:));% 
                    vel(k,:)=iwt(i).*vel(k,:)+aa(k,:); 
                    vel(k,:)=(vel(k,:)>mv).*mv+(vel(k,:)<=mv).*vel(k,:); 
                    vel(k,:)=(vel(k,:)<(-mv)).*(-mv)+(vel(k,:)>=(-mv)).*vel(k,:);
                    pos(k,:)=pos(k,:)+vel(k,:); 

                    if (sum(pos(k,:)>VRmax(k,:))+sum(pos(k,:)<VRmin(k,:)))==0;% if the position of the updated particle is not out of scope, gbest and pbest are updated.
                        %e(k)=feval(fhd,pos(k,:),varargin{:});
                        e(k)=fobj(pos(k,:));
                        fitcount=fitcount+1;
                        tmp=(pbestval(k)<=e(k));
                        if tmp==1
                            stay_num(k)=stay_num(k)+1;
                        end
                        temp=repmat(tmp,1,D);
                        pbest(k,:)=temp.*pbest(k,:)+(1-temp).*pos(k,:);
                        pbestval(k)=tmp.*pbestval(k)+(1-tmp).*e(k);%update the pbest
                        if pbestval(k)<gbestval
                            gbest=pbest(k,:);
                            gbestval=pbestval(k);
                            gbestrep=repmat(gbest,ps,1);%update the gbest
                        end
                    end   % if the position of the updated particle is out of scope, gbest and pbest are not updated.
                        
            end % corresponding 64 line

% if round(i/100)==i/100
%     plot(pos(:,D-1),pos(:,D),'b*');hold on;
%     for k=1:floor(D/2)
%         plot(gbest(:,2*k-1),gbest(:,2*k),'r*');
%     end
%     hold off
%     title(['PSO: ',num2str(i),' generations, Gbestval=',num2str(gbestval)]);  
%     axis([VRmin(1,D-1),VRmax(1,D-1),VRmin(1,D),VRmax(1,D)])
%     drawnow
% end
allgbestval(i)=gbestval;

        if fitcount>=Max_FES
            break;
        end
        if mod(i,100)==0
            disp(['iter: ',num2str(i),'   gebestVal: ',num2str(gbestval)]);
        end
        if (i==me)&(fitcount<Max_FES)% when the position is out of the scope, the fitness is not assessed. So, fitcout may be less than me.
            i=i-1;
        end

end %corresponding 62 line
gbestval;
gbest;
allgbestval;
fitcount;
end


