%%
%**************************************************************************
%*******************ע�⣬�˳�������Ӧ��ֵ��С��Ϊ����***********************
%**************************************************************************
clc
clear all
format long
%% ������ʼ��
for cf=1:10
X_min=0;                                                     %�����С
X_max=3;

P_ini=10;                                                   %��ʼ�Ӳݸ���
P_max=20;                                                   %����Ӳݸ���

D=20;                                                        %�����ά��
iter_max= 200;                                                %��������

stepLength_ini=5;                                            %��ʼ��׼��
stepLength_final=1;                                          %���ձ�׼��
                        
seed_max=3;                                                  %���������
seed_min=0;                                                  %��С������
n=3;                                                         %����������

L = 11;
penalty = 10^(30);

                                             
tic
%% ��ʼ��
X_ini=X_min+round((X_max-X_min).*rand(P_ini,D));                %Dά�ռ����漴�ֲ���G_SIZE�����н�
X=X_ini;
f=fit(X);                                                   %��ȡ��Ӧ��
rank=f(:,5);
X=[X,rank];
X=sortrows(X,D+1);
%% ����
for iter = 1:iter_max
    N=size(X,1);
    BestFitness = X(1,D+1);
    WorstFitness = X(end,D+1);
    avg_fit=sum(X(:,D+1))/N;
   %%  �������Ӹ�������������
    stepLength_now=(iter_max-iter)^n*(stepLength_ini-stepLength_final)/(iter_max)^n+stepLength_final;   % ���㲽��������ǰ��׼��
    num=(seed_max-seed_min)*(X(:,D+1)-WorstFitness)/(BestFitness-WorstFitness)+seed_min;        % ��������Ӳ������������Ӹ���
    num=floor(num);    %����ȡ��
   
 
    X1=[];                                                    %ȫ������
    for i=1:N
        for j = 1:num(i)
            weed_new = mod(X(i,1:D)+round(stepLength_now*trnd(1,1,D)),4);
            X1=[X1;weed_new];                                 %�������������Ӵ�����X1��
        end
    end
    
    X2=[X(:,1:D);X1];
    N1=size(X2,1);
    
    
    f1=fit(X2);                                                   %��ȡ��Ӧ��
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
        f=fit(X);                                                   %��ȡ��Ӧ��
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
    f=fit(X);                                                   %��ȡ��Ӧ��
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
toc  %��ʾ��������ʱ�� 

