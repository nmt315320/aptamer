%%
%**************************************************************************
%*******************注意，此程序以适应度值最小作为最优***********************
%**************************************************************************
clc
clear all
format long
%% 参数初始化
for cf=1:10
X_min=0;                                                     %区间大小
X_max=3;

P_ini=10;                                                   %初始杂草个数
P_max=20;                                                   %最大杂草个数

D=20;                                                        %问题的维数
iter_max= 200;                                                %迭代次数

stepLength_ini=5;                                            %初始标准差
stepLength_final=1;                                          %最终标准差
                        
seed_max=3;                                                  %最大种子数
seed_min=0;                                                  %最小种子数
n=3;                                                         %非线性因子

L = 11;
penalty = 10^(30);

                                             
tic
%% 初始化
X_ini=X_min+round((X_max-X_min).*rand(P_ini,D));                %D维空间中随即分布的G_SIZE个可行解
X=X_ini;
f=fit(X);                                                   %求取适应度
rank=f(:,5);
X=[X,rank];
X=sortrows(X,D+1);
%% 进化
for iter = 1:iter_max
    N=size(X,1);
    BestFitness = X(1,D+1);
    WorstFitness = X(end,D+1);
    avg_fit=sum(X(:,D+1))/N;
   %%  计算种子个数并产生种子
    stepLength_now=(iter_max-iter)^n*(stepLength_ini-stepLength_final)/(iter_max)^n+stepLength_final;   % 计算步长，即当前标准差
    num=(seed_max-seed_min)*(X(:,D+1)-WorstFitness)/(BestFitness-WorstFitness)+seed_min;        % 计算各个杂草所产生的种子个数
    num=floor(num);    %向下取整
   
 
    X1=[];                                                    %全部种子
    for i=1:N
        for j = 1:num(i)
            weed_new = mod(X(i,1:D)+round(stepLength_now*trnd(1,1,D)),4);
            X1=[X1;weed_new];                                 %将产生的所有子代存在X1中
        end
    end
    
    X2=[X(:,1:D);X1];
    N1=size(X2,1);
    
    
    f1=fit(X2);                                                   %求取适应度
    rank1=f1(:,5);
    X2=[X2,rank1];
    X2=sortrows(X2,D+1);

    X=[];  
    if N1>P_max
        X2 = Niche(N1,D,L,penalty,X2);
        [X2,q] = sortrows(X2,D+1);
        X=X2(1:P_max,1:D);
        f=[];
        rank=[];
        f=fit(X);                                                   %求取适应度
        rank=f(:,5);
        X=[X,rank];
        X=sortrows(X,D+1);
        X=X(:,1:D);
    else
        X=X2(:,1:D);
    end
    N2=size(X,1);
    
    
   

    X=repair(N2,D,X);
    f=[];
    rank=[];
    f=fit(X);                                                   %求取适应度
    rank=f(:,5);
    X=[X,rank];
    X=sortrows(X,D+1);
end

for i=1:7
    X3(i,:)=X(i,1:D);
end
f3=fit(X3);
D=DNAcode2(X3);
[f3(:,5),f3(:,6)] = GCTmBioBox(D);
R=[X3,f3];
eval(['save(''','Result',num2str(cf),'.mat''',',','''R''',');']);
end
toc  %显示程序运行时间 

