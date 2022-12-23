%---------------------����ϸ����ʳ�㷨��DNA����������------------------------------%
%-----------------------���ߣ�ҦҢ     ��λ����������ѧ--------------------------------%

%----------------------------------ģ��ϸ����ʳ����--------------------------------------%
%----------------------------------ģ��ϸ����ʳ����--------------------------------------%
%----------------------------------ģ��ϸ����ʳ����--------------------------------------%
%----------------------------------ģ��ϸ����ʳ����--------------------------------------%
%----------------------------------ģ��ϸ����ʳ����--------------------------------------%
close all;
clc;
%----------------------
C=0.01;    %ͻ�����

%Ne=20;
%Np=20;
%Ns=10;��ʼ��------------------------%
% ------- ��ʼ�� ----------%
Nre=10;    %���ƴ���
Nc=10;    %��������������������
Ped=0.1;  %�������
%D=5;

BFOA
evaluate(DNA)




%-----------------------����------------------------%
function evaluate(dna)
%��Լ������һ����Ȩ��������E
%Hanming=hanming(DNA(i,:))
%Similarity=similarity(DNA(i,:))
%H_measure=h_measure(DNA(i,:))
%[Similarity1,H_measure1]=HmSm(DNA(i,:))
%Continuity=continuity(DNA(i,:))
%GC=gccontent(DNA(i,:))
%Tm=TmBioBox(DNA(i,:))
%Hairpinnum= Hairpin(DNA(i,:))
%�Ʒֹ�������
%7�����������У�Tm40�֣�hanmingռ10�֡�simռ15�֡�H_mռ15�֡�conռ20�֣�����100��
%1��GC������50%��hairpin������0���ܷ�ֱ�Ӽ�Ϊ0��
%2��Tm����<=60��=40�֣���60~70��40���𼶵ݼ�����������ɢ�ݼ�������>=70��=0��
%3��hanming����=20��=10�֣���20~16��ÿλ��2�֣���<15��=0��
%4��sim������20����0~10��ÿλ��1.5,���ԣ���>10��=0��
%5��H_m������20����0~10��ÿλ��1.5,���ԣ���>10��=0��
%6��con����=0��=20�֣���=1��=15�֣���=2��=5�֣���>=3��=0��
 goal=double(20);
 for i=1:20
     goal(i)=0;
 end
 for i=1:20
     
     if similarity(dna(i,:),dna(:,:))/20<=10
         goal(i)=goal(i)+(10-similarity(dna(i,:),dna(:,:))/20)*1.5;
     else
         goal(i)=goal(i)+0;
     end
     if h_measure(dna(i,:),dna(:,:))/20<=10
         goal(i)=goal(i)+(10-h_measure(dna(i,:),dna(:,:))/20)*1.5;
     else
         goal(i)=goal(i)+0;
     end
     if continuity(dna(i,:))==0
         goal(i)=goal(i)+20;
     elseif continuity(dna(i,:))==1
         goal(i)=goal(i)+15;
     elseif continuity(dna(i,:))==2
         goal(i)=goal(i)+5;
     else 
         goal(i)=goal(i)+0;
     end
     if TmBioBox(dna(i,:))<=60
         goal(i)=goal(i)+40;
     elseif  TmBioBox(dna(i,:))>60 && TmBioBox(dna(i,:))<=70
         goal(i)=goal(i)+((70-TmBioBox(dna(i,:)))/10)*40;
     else
         goal(i)=goal(i)+0;
     end
     if gccontent(dna(i,:))==0.5
         goal(i)=goal(i);
     else
         goal(i)=0;
     end
     if Hairpin(dna(i,:))==0
          goal(i)=goal(i);
     else
         goal(i)=0;
     end
     
     fprintf ('%d��DNA�ĵ÷�=%.2f\n',i,goal(i))
 end   
end


function h=hanming(x,y)
h=0;
for i=1:20
if x(i)==y(i)
    h=h+1;
end
end
h=20-h;
end
%-----------------------------------------------------------------------------------------%

%--------------------------------------���ƶ�Լ��---------------------------------------%
function S=similarity(x,y)
s=double(38);
for i=1:38
    s(i)=0;
end
for k=1:19  %k��ʾy��'��'�ƶ�����
    y=circshift(y,k);
    for i=k:19  
        if x(i)==y(i)
            s(k)=s(k)+1;
        end
    end
end
for kk=1:19 %kk��ʾy��'��'�ƶ�����(��x�����ƶ����y�����ƶ�)
    x=circshift(x,kk);
    for i=kk:19  
        if x(i)==y(i)
            s(19+kk)=s(kk)+1;
        end
    end
end
[S,sp]=max(s); %sp��ʾ����������Ƶ�λ�ã�19֮ǰΪy���ƣ�20֮��Ϊy���ƣ�
end
%-----------------------------------------------------------------------------------------%

%--------------------------------------H-measureԼ��---------------------------------------%
function hm=h_measure(x,y)
s=double(38);
for i=1:38
    s(i)=0;
end
for yc=1:20   %��yȡ��
    switch y(yc)
        case 'A'
            y(yc)='T';
        case 'T'
            y(yc)='A';
        case 'G'
            y(yc)='C';
        case 'C'
            y(yc)='G';
    end
end
for k=1:19  %k��ʾy��'��'�ƶ�����
    y=circshift(y,k);
    for i=k:19  
        if x(i)==y(i)
            s(k)=s(k)+1;
        end
    end
end
for kk=1:19 %kk��ʾy��'��'�ƶ�����(��x�����ƶ����y�����ƶ�)
    x=circshift(x,kk);
    for i=kk:19  
        if x(i)==y(i)
            s(19+kk)=s(kk)+1;
        end
    end
end
[hm,hmp]=max(s);%hmp��ʾ����������Ƶ�λ�ã�19֮ǰΪy���ƣ�20֮��Ϊy���ƣ�
end
%-----------------------------------------------------------------------------------------%

%-----------------------------------�ظ��ԣ������ԣ�Լ��--------------------------------%
function con=continuity(x)
c1=1;
c2=1;
cmax=1;
x(21)=0;  %��֤ĩβ��ͬ
for i=1:20
if x(i)==x(i+1) 
    c1=c1+1;
else
    c2=c1;
    c1=1;
if c2>=cmax 
    cmax=c2;
end
end
end
con=cmax;
end
%----------------------------------------------------------------------------------%

%-----------------------------------GC����Լ��------------------------------------%
function gc=gccontent(x)
p=0;
for i=1:20
    if x(i)=='G' || x(i)=='C'
        p=p+1;
    end
end
p=(p/20)*100;
gc=num2str(p)+'%';
end
%----------------------------------------------------------------------------------%

%-----------------------------------------TmԼ��----------------------------------%
function [Tm,GC] = TmBioBox(x)
SaltValue=1; % Shin2005 ��Ũ��1M  Bioinformatics Toolbox Ĭ�� 0.05
PrimerconcValue=10^(-8); %DNA����Ũ��
m=size(x,1);% DNA���м������
GC=zeros(m,1);
Tm=GC;
for i=1:m 
    oligoprops=oligoprop(x(i,:),'Salt', SaltValue,'Primerconc', PrimerconcValue);
   Tm(i)=round(oligoprops.Tm(5)*10000)/10000; %% Tmֵ������С�������λ
end
end
%----------------------------------------------------------------------------------------%

%------------------------------------------����Լ��---------------------------------------%
function Values = Hairpin( DNAs )
% ����DNAs(�����������)�ķ�����������Values������
% function IntValue = Hairpin( IntDNA, R_min,P_min )
%��Shin2005,����IntDNA�ķ�����
R_MIN=6;
P_MIN=6;
PINLEN=3;
if nargin==1
    R_min=R_MIN;
    P_min=P_MIN;
end
[m,l]=size(DNAs);
Values=zeros(m,1);
for k=1:m
    IntDNA=DNAs(k,:);
    for p=P_min:(l-R_min)/2
        for r=R_min:(l-2*p)
            for i=0:(l-2*p-r)
                Sigma_bp=0;
                for j=0:pinlen(p,r,i,IntDNA)-1
                    Sigma_bp=Sigma_bp+bp(IntDNA(p+i-j),IntDNA(p+i+r+1+j));
                end
                Values(k)=Values(k)+T(Sigma_bp,pinlen(p,r,i,IntDNA)/ PINLEN);
            end
        end
    end
end
end

function minValue = pinlen( p,r,i,IntDNA )
%��������������min(p+i,l-(p+i+r))
l=size(IntDNA,2);
minValue=min(p+i,l-p-i-r);
end 

function BoolValue = bp( IntDNAa,IntDNAb )
%intDNAa��intDNAb�Ƿ񻥲���Y����1,N����0;IN����������DNA����{ACGT}->{2031};OUT������ֵ{0,1}
l=size(IntDNAa,2);
m=size(IntDNAb,2);
if l==m
    temp=3-IntDNAa;
    if sum(temp==IntDNAb)==l
        BoolValue=1;%'����'
    else
        BoolValue=0;%'������'
    end
else
    BoolValue=0;
    error('DNA���Ȳ�ͬ');
end
end

function IntValue= T(i,j)
%i>j����i�����򷵻�0
if i>j 
    IntValue=i;
else
    IntValue=0;
end
end