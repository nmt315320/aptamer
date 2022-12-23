%---------------------基于细菌觅食算法的DNA计算编码设计--------------------------%
%-----------------------作者：姚尧     单位：大连理工大学----------------------------%
close all;
clc;
% ------- 初始化 ----------%
Nr=20;    %复制次数
Nc=20;    %迭代次数（趋化次数）
Ped=0.9;  %消灭概率
C=0.01;    %突变概率

%Ne=20;
%Np=20;
%Ns=10;
%D=5;

%-----随机生成20个碱基DNA序列作为细菌-----%
%-----定义碱基---0=A；1=G；2=T；3=C；---%
dna=randi([0,3],2,20);
DNAx=char(20);
DNAy=char(20);
for i=1:2
    for ii=1:20
        if i==1         %DNAx
switch dna(1,ii)
    case 0
        DNAx(ii)='A';
        num2str(DNAx(ii));
    case 1
        DNAx(ii)='G';
        num2str(DNAx(ii));
    case 2
        DNAx(ii)='T';
        num2str(DNAx(ii));
    case 3
        DNAx(ii)='C';
        num2str(DNAx(ii));
end
        end
         if i==2     %DNAy
switch dna(2,ii)
    case 0
        DNAy(ii)='A';
        num2str(DNAy(ii));
    case 1
        DNAy(ii)='G';
        num2str(DNAy(ii));
    case 2
        DNAy(ii)='T';
        num2str(DNAy(ii));
    case 3
        DNAy(ii)='C';
        num2str(DNAy(ii));
end
        end
    end
end
DNA=[DNAx',DNAy']'
%-----------------------------------------------%

%--------------约束检测----------------%
Hanming=hanming(DNAx,DNAy)
Similarity=similarity(DNAx,DNAy)
H_measure=h_measure(DNAx,DNAy)
Continuity=continuity(DNAx)
GC=gccontent(DNAx)
Tm=TmBioBox(DNAx)
Hairpinnum= Hairpin(DNAx)
%---------------------------------------%




%--------------------------------------各种约束-----------------------------------------%
%------------------------------------汉明距离约束--------------------------------------%
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

%--------------------------------------相似度约束---------------------------------------%
function S=similarity(x,y)
s=double(38);
for i=1:38
    s(i)=0;
end
for k=1:19  %k表示y向'右'移动距离
    y=circshift(y,k);
    for i=k:19  
        if x(i)==y(i)
            s(k)=s(k)+1;
        end
    end
end
for kk=1:19 %kk表示y向'左'移动距离(用x向右移动替代y向左移动)
    x=circshift(x,kk);
    for i=kk:19  
        if x(i)==y(i)
            s(19+kk)=s(kk)+1;
        end
    end
end
[S,sp]=max(s); %sp表示出现最大相似的位置（19之前为y右移，20之后为y左移）
end
%-----------------------------------------------------------------------------------------%

%--------------------------------------H-measure约束---------------------------------------%
function hm=h_measure(x,y)
s=double(38);
for i=1:38
    s(i)=0;
end
for yc=1:20   %对y取补
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
for k=1:19  %k表示y向'右'移动距离
    y=circshift(y,k);
    for i=k:19  
        if x(i)==y(i)
            s(k)=s(k)+1;
        end
    end
end
for kk=1:19 %kk表示y向'左'移动距离(用x向右移动替代y向左移动)
    x=circshift(x,kk);
    for i=kk:19  
        if x(i)==y(i)
            s(19+kk)=s(kk)+1;
        end
    end
end
[hm,hmp]=max(s);%hmp表示出现最大相似的位置（19之前为y右移，20之后为y左移）
end
%-----------------------------------------------------------------------------------------%

%-----------------------------------重复性（连续性）约束--------------------------------%
function con=continuity(x)
c1=1;
c2=1;
cmax=1;
x(21)=0;  %保证末尾不同
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

%-----------------------------------GC含量约束------------------------------------%
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

%-----------------------------------------Tm约束----------------------------------%
function [Tm,GC] = TmBioBox(x)
SaltValue=1; % Shin2005 盐浓度1M  Bioinformatics Toolbox 默认 0.05
PrimerconcValue=10^(-8); %DNA分子浓度
m=size(x,1);% DNA序列碱基数量
GC=zeros(m,1);
Tm=GC;
for i=1:m 
    oligoprops=oligoprop(x(i,:),'Salt', SaltValue,'Primerconc', PrimerconcValue);
   Tm(i)=round(oligoprops.Tm(5)*10000)/10000; %% Tm值，保留小数点后四位
end
end
%----------------------------------------------------------------------------------------%

%------------------------------------------发夹约束---------------------------------------%
function Values = Hairpin( DNAs )
% 计算DNAs(整数编码矩阵)的发卡数，返回Values列向量
% function IntValue = Hairpin( IntDNA, R_min,P_min )
%依Shin2005,计算IntDNA的发卡数
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
%辅助函数，返回min(p+i,l-(p+i+r))
l=size(IntDNA,2);
minValue=min(p+i,l-p-i-r);
end 

function BoolValue = bp( IntDNAa,IntDNAb )
%intDNAa与intDNAb是否互补，Y返回1,N返回0;IN：整数编码DNA向量{ACGT}->{2031};OUT：整数值{0,1}
l=size(IntDNAa,2);
m=size(IntDNAb,2);
if l==m
    temp=3-IntDNAa;
    if sum(temp==IntDNAb)==l
        BoolValue=1;%'互补'
    else
        BoolValue=0;%'不互补'
    end
else
    BoolValue=0;
    error('DNA长度不同');
end
end

function IntValue= T(i,j)
%i>j返回i，否则返回0
if i>j 
    IntValue=i;
else
    IntValue=0;
end
end
%------------------------------------------------------------------------------------------%