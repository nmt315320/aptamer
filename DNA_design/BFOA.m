%---------------------����ϸ����ʳ�㷨��DNA����������--------------------------%
%-----------------------���ߣ�ҦҢ     ��λ������������ѧ----------------------------%

%--------------------------------��ǰ������DNA���ɼ�Լ��---------------------------%
%--------------------------------��ǰ������DNA���ɼ�Լ��---------------------------%
%--------------------------------��ǰ������DNA���ɼ�Լ��---------------------------%
%--------------------------------��ǰ������DNA���ɼ�Լ��---------------------------%
%--------------------------------��ǰ������DNA���ɼ�Լ��---------------------------%
%ʹ��˵����1.����BFOA��ʼ�������������ơ����������ʡ����������ȣ�                            %
%                2.������У��ȴ���������ʾ������ʱ�䣬��Ϊ���γ��������                          %
%                   �мǣ�����������ֻ��������ʾ��Ž�չ�ģ����Ը��ݽ�������Ϣȥ�ȵ�裻    %
%                   �мǣ�������������ʱ��󣬹رս��������ڣ�                                             %
%                3.������ʱ�����������ΪnumSave�����У��������Ƶ�finalchoose�ӳ����У� %
%                  ����˵����finalchoose�У���ת���鿴                                                         %


close all;
clc;
tic
%----------------------��ʼ��------------------------%
% ------- ��ʼ�� ----------%
Nre=9;    %���ƴ���
Nc=10;    %��������������������
Ped=0.1;  %�������
numDNA=20; %����DNA����
numSave=10;%���ձ���DNA����
%Ne=20;
%Np=20;
%Ns=10;
%D=5;
diary('E:\��ʿ�������\������\DNA_design\record.txt');
delete('E:\��ʿ�������\������\DNA_design\record.txt');
diary('E:\��ʿ�������\������\DNA_design\record.txt');

%-----�������20�����DNA������Ϊϸ��-----%
%-----������---0=A��1=G��2=T��3=C��---%
Dna=randi([0,3],numDNA,20);
DNA=char(numDNA,20);
%d1=['GTACTCTACTACAGCTGCAG';'GCACGTAGATGCAAAAGTCG';'ACTAGTGTGTATACGAGTGC'];
%evaluate(d1,numDNA,numSave)

for i=1:numDNA
    for ii=1:20
        switch Dna(i,ii)
            case 0
                DNA(i,ii)='A';
                num2str(DNA(i,ii));
            case 1
                DNA(i,ii)='G';
                num2str(DNA(i,ii));
            case 2
                DNA(i,ii)='T';
                num2str(DNA(i,ii));
            case 3
                DNA(i,ii)='C';
                num2str(DNA(i,ii));
        end
        
    end
end
fprintf('���ɵ�DNA����Ϊ��\n')
for i=1:numDNA
    if i<10
        fprintf('  %d:::     5''-%s-3''\n',i,DNA(i,:))
    else
        fprintf('%d:::     5''-%s-3''\n',i,DNA(i,:))
    end
end

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-%
%***********************************************%
%***********************************************%
% DNA=['GTACTCTACTACAGCTGCAG';'GCACGTAGATGCAAAAGTCG';'ACTAGTGTGTATACGAGTGC';'GAGCAGTGTACTAGACTGTG';'GTACACACAGCAGACATCTC';
% 'GCTCTGTAGTGACGAGTATC';'GCATATGTCTCTCGCAGATG';'ATCGAGTCTGAGAGACACTG';'GATGCGAGCTACTAGTGATG';'TAGACGACTACGCACTACTC';
% 'CACGAGTATAGTCAGCACTG';'CCTATTACGTGGACTATCCG';'ACGCAGTGCTATACTAGCAG';'ACGAGTGTGCGTGTATGTAG';'TGACGCTCTATAGTGCACTG';
% 'AGAGAGACGTCGACTCAGTA';'CTATATACGAGAGCGCACGA';'TGATGTAGATGCTGAGCTCG';'GAGCAGTGTACTAGACTGTG';'GTACACACAGCAGACATCTC';
% 'GCTCTGTAGTGACGAGTATC';'GCATATGTCTCTCGCAGATG';'ATCGAGTCTGAGAGACACTG';'GATGCGAGCTACTAGTGATG';'TAGACGACTACGCACTACTC']

e=0;
w=waitbar(0,'��ʼ'); 
for i=1:Nre    %���ƴ���
    if i<Nre
    process=i/Nre;    
    waitbar(i/Nre,w,['ū���Ѿ�����' num2str(100*process) '%�ˣ��������ٺȵ��']);
    elseif i==Nre
    waitbar(1,w,'���˳��ˣ����Ҫ���ˣ�����');   
    end
    diary on;
    analysis(DNA,numDNA);
    evaluate(DNA,numDNA,numSave);
    if i>1
        fprintf('����ɵ�%d�θ���--------------------------------------����ɵ�%d�θ���\n',i-1,i-1)
        fprintf('�����ĵȴ�����ʣ%d�θ���\n',Nre-i+1)
    end
    DNA=Chemotaxis(DNA,Nc,numDNA);
    fprintf('�޸ĺ�--------------------------------------�޸ĺ�\n')
    fprintf('�޸ĺ�--------------------------------------�޸ĺ�\n')
    fprintf('��%d���޸ĺ�DNA����Ϊ��\n',i)
    for iii=1:numDNA
        if iii<10
            fprintf('  %d:::     5''-%s-3''\n',iii,DNA(iii,:))
        else
            fprintf('%d:::     5''-%s-3''\n',iii,DNA(iii,:))
        end
    end
    analysis(DNA,numDNA);
    BetterDNAindex=evaluate(DNA,numDNA,numSave);
    if i<Nre %���1�β�����
        [DNA,ped]=reproduce(DNA(BetterDNAindex(1:numDNA/2),:));%�÷ָߵ�һ��DNA�����ƣ������̭
    end
    if ped<Ped && i<Nre %���һ�β�����
        e=e+1;
        DNA=eliminate(DNA);
        fprintf('��������%d��--------------------------------------�������\n',e)
        fprintf('��������%d��--------------------------------------�������\n',e)
    end
    fprintf('�Ѹ���%d��--------------------------------------�Ѹ���%d��\n',i,i)
    fprintf('�Ѹ���%d��--------------------------------------�Ѹ���%d��\n',i,i)
end
fprintf('���һ�θ��������--------------------------------------���һ�θ��������\n')
fprintf('*************************�����������������DNA����************************\n')
BetterDna=evaluate(DNA,numDNA,numSave);
for ii=1:numSave
    if ii<10
        fprintf('  %d��%d��:::     5''- %s -3''\n',ii,BetterDna(ii),DNA(BetterDna(ii),:))
    else
        fprintf('%d��%d��:::     5''- %s -3''\n',ii,BetterDna(ii),DNA(BetterDna(ii),:))
    end
    
end
diary off;
toc
%disp(['����ʱ��: ',num2str(toc)]);
%***********************************************%
%***********************************************%
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-%

%--------------Լ�����----------------%
function analysis(dna,numDNA)
for i=1:numDNA
    fprintf('��%d��DNA��Լ�������',i)
    fprintf(dna(i,:))
    fprintf('\n')
    %fprintf('Hanming=%d\n',hanming(dna(i,:),dna(:,:)))
    %s=dna(:,:);
    %s(i,:)=[];
    [sim,~]=similarity(dna(i,:),dna(:,:),numDNA,i);
    [hme,~]=h_measure(dna(i,:),dna(:,:),numDNA);
    fprintf('Similarity=%d    \n',sim)
    fprintf('H_measure=%d      \n',hme)
    fprintf('Continuity=%d\n',continuity(dna(i,:),numDNA))
    %fprintf('Continuity111=%d\n',Continuity(DNA(i,:)))
    fprintf('GC_content=%.2f\n',gccontent(dna(i,:),numDNA))
    fprintf('Tm=%.2f\n',TmBioBox(dna(i,:)))
    fprintf('Hairpinnum=%d\n',Hairpin(dna(i,:)))
    %Hanming=hanming(DNA(i,:))
    %Similarity=similarity(DNA(i,:))
    %H_measure=h_measure(DNA(i,:))
    %[Similarity1,H_measure1]=HmSm(DNA(i,:))
    %Continuity=continuity(DNA(i,:))
    %GC=gccontent(DNA(i,:))
    %Tm=TmBioBox(DNA(i,:))
    %Hairpinnum= Hairpin(DNA(i,:))
end
end
%---------------------------------------%
function fianl(dna,numDNA,numSave,fianlbest)
for i=1:numSave
    fprintf('��'+i+'('+fianlbest(i)+')'+'��DNA��Լ�������')
    fprintf(dna(fianlbest(i),:))
    fprintf('\n')
    %s=dna;
    %s(i,:)=[];
    [sim,~]=similarity(dna(fianlbest(i),:),dna(:,:),numDNA,i);
    [hme,~]=h_measure(dna(fianlbest(i),:),dna(:,:),numDNA);
    fprintf('Similarity=%d      \n',sim)
    fprintf('H_measure=%d      \n',hme)
    fprintf('Continuity=%d\n',continuity(dna(fianlbest(i),:),numDNA))
    fprintf('GC_content=%.2f\n',gccontent(dna(fianlbest(i),:),numDNA))
    fprintf('Tm=%.2f\n',TmBioBox(dna(fianlbest(i),:)))
    fprintf('Hairpinnum=%d\n',Hairpin(dna(fianlbest(i),:)))
end
end


%--------------------------------------����Լ��-----------------------------------------%
%------------------------------------��������Լ��--------------------------------------%
function h=hanming(x,y)
han=char(20);
for i=1:20
    han(i)=0;
end
for i=1:20
    for ii=1:20
        for iii=1:20
            if x(i,ii)~=y(i,iii)
                han(i)=han(i)+1;
            end
        end
    end
end

h=max(han);
end
%-----------------------------------------------------------------------------------------%

%--------------------------------------���ƶ�Լ��---------------------------------------%
function [S,Smax]=similarity(xo,yo,numDNA,flag)
s=[];
for i=1:40
    s(i)=0;
end
smax=[];
for i=1:numDNA
    smax(i)=0;
end

stotal=0;
Smax=0;
for ii=1:numDNA
    if flag~=ii
        for k=1:19  %k��ʾy��'��'�ƶ�����
            x=xo;
            y=circshift(yo(ii,:),k);
            for i=k:20
                if  x(i)==y(i)
                    s(k)=s(k)+1;
                end
            end
        end
        for kk=1:19 %kk��ʾy��'��'�ƶ�����(��x�����ƶ����y�����ƶ�)
            y=yo;
            x=circshift(xo,kk);
            for i=kk:20
                if  x(i)==y(ii,i)
                    s(19+kk)=s(kk+19)+1;
                end
            end
        end
        for i=1:20
            if  xo(i)==yo(ii,i)
                s(40)=s(40)+1;
            end
        end
        [ss,~]=max(s);%~��ʾ����������Ƶ�λ�ã�20֮ǰΪy���ƣ�21֮��Ϊy���ƣ�
        s(:,:)=0;
        stotal=stotal+ss;
        smax(ii)=ss;
    end
end
[Smax,~]=max(smax);
S=stotal;
end
%-----------------------------------------------------------------------------------------%

%--------------------------------------H-measureԼ��---------------------------------------%
function [hm,Hmax]=h_measure(xo,yo,numDNA)
s=[];
for i=1:40
    s(i)=0;
end
hmax=[];
for i=1:numDNA
    hmax(i)=0;
end
htotal=0;
for ci=1:numDNA   %��yȡ��
    for cii=1:20
    switch yo(ci,cii)
        case 'A'
            yo(ci,cii)='T';
        case 'T'
            yo(ci,cii)='A';
        case 'G'
            yo(ci,cii)='C';
        case 'C'
            yo(ci,cii)='G';
    end
    end
end
for ii=1:numDNA
    for k=1:19  %k��ʾy��'��'�ƶ�����
        x=xo;
        y=circshift(yo(ii,:),k);
        for i=k:20
            if  x(i)==y(i)
                s(k)=s(k)+1;
            end
        end
    end
    for kk=1:19 %kk��ʾy��'��'�ƶ�����(��x�����ƶ����y�����ƶ�)
        y=yo;
        x=circshift(xo,kk);
        for i=kk:20
            if  x(i)==y(ii,i)
                s(19+kk)=s(kk+19)+1;
            end
        end
    end
    for i=1:20
        if  xo(i)==yo(ii,i)
            s(40)=s(40)+1;
        end
    end
    [ss,~]=max(s);%hmp��ʾ����������Ƶ�λ�ã�20֮ǰΪy���ƣ�21֮��Ϊy���ƣ�
    s(:,:)=0;
    htotal=htotal+ss;
    hmax(ii)=ss;
end
[Hmax,~]=max(hmax);
hm=htotal;
end

%-----------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------%
function [Hm,Sim]=HmSm(DNAs)
% �����������ۣ�����DNA����DNA1�������DNAs�����еĵ�Similarity��H-measure��������Shin2005������ͬ������2������������DNA1��DNAs�е�ĳ��������ͬ������SimilarityΪ0
% DNA1�ĳ�����DNAs��������ͬ,��ͬ���е�Similarity=0
% Similarity����Shin2005��H-measure����Shin2005,��c0720_Shin2005TSP_MEA��
% global Gap;
[m,l]=size(DNAs);
Gap1=round(l/4);  % ȡ���е�1/4ΪGap����Shin2005���ݻ�����ͬ
for p=1:m
    SimValue=zeros(m,1);
    HmValue=zeros(m,1);
    DNA1=DNAs(p,:);
    ReverseDNA1=seqreverse(DNA1);%IntDNAxȡ��
    for j=1:m
        if sum(DNA1==DNAs(j,:))==l  % ��ͬ��DNA���в��Ƚ�Similarity,���Ƚ�H-measure
            for g=0:Gap1  %  ȡ���е�1/4ΪGap
                tempIntDNAy=[DNAs(j,:), zeros(1,g)+5, DNAs(j,:)];
                for i=-l+1:l-1
                    ShiftValue=shift(tempIntDNAy,i);
                    currentHm=h_dis(ReverseDNA1,ShiftValue)+h_con(ReverseDNA1,ShiftValue);
                    HmValue(j)=max(HmValue(j),currentHm);
                end
            end
        else                        % ͬʱ�Ƚ�Similarity��H-measure
            for g=0:Gap1  %  ȡ���е�1/4ΪGap
                tempIntDNAy=[DNAs(j,:), zeros(1,g)+5, DNAs(j,:)];
                for i=-l+1:l-1
                    ShiftValue=shift(tempIntDNAy,i);
                    currentSim=s_dis(DNA1,ShiftValue)+s_con(DNA1,ShiftValue);
                    currentHm=h_dis(ReverseDNA1,ShiftValue)+h_con(ReverseDNA1,ShiftValue);
                    SimValue(j)=max(SimValue(j),currentSim);
                    HmValue(j)=max(HmValue(j),currentHm);
                end
            end
        end
    end
    Hm(p)=sum(HmValue);
    Sim(p)=sum(SimValue);
end
end

function IntDNAShifted= shift( IntDNA_nb,i )
%��IntDNA_nbƽ��iλ��i>0���ƣ�i<0,����;��λ��'-'(5)����
if i==0
    IntDNAShifted=IntDNA_nb;
    return;
end
l=size(IntDNA_nb,2);%DNA����
temp=zeros(1,abs(i))+5;
if (0<i) && (i<l)   %����iλ����i��'-'
    IntDNAShifted=[temp, IntDNA_nb(1:l-i)];
    return;
end
if (i<0) && (i>-l)  %����iλ����i��'-'
    IntDNAShifted=[IntDNA_nb(abs(i)+1:l),temp];
    return;
end
if abs(i)>=l %ȫ���滻Ϊ'-'
    IntDNAShifted=zeros(1,l)+5;
    return;
end
end

function IntValue= h_dis( IntDNAx,IntDNAy,H_dis)
%H_measure������һ�H_disΪ�û�ָ������;H_disΪ0..1
global H_DIS;
H_DIS = 0.17;
if nargin==2
    H_dis=H_DIS;
end
Sigma_bp=0;
l=size(IntDNAx,2);
for i=1:l
    Sigma_bp=Sigma_bp+bp(IntDNAx(i),IntDNAy(i));
end
temp=H_dis*length_nb(IntDNAy)/2;
IntValue=T(Sigma_bp,temp);
end

function LenNoBlank= length_nb(IntDNA_nb)
%������ĸ��{A C G T -}->{1 2 3 4 5}���ַ����з�-���ַ���
LenNoBlank=sum(IntDNA_nb~=5);
end

function IntValue = h_con( IntDNAx,IntDNAy,H_con )
%H_measure�����ڶ��H_conΪ�û�ָ������;H_conΪ1..Length(DNA)
global H_CON;
H_CON=6;
if nargin==2
    H_con=H_CON;
end
IntValue=0;
l=size(IntDNAx,2);
for i=1:l
    IntValue=IntValue+T(cbp(IntDNAx,IntDNAy,i),H_con);
end
end

function IntValue = cbp(IntDNAx,IntDNAy,i)
%�ӵ�i��λ�ÿ�ʼ����Ϊc��IntDNAx��IntDNAy�Ĺ��������Ӵ��ĳ���;IN:IntDNAx��IntDNAy,��������DNA,i,��ʼλ��;OUT:IntValue,������ִ�����
IntValue=0;
l=size(IntDNAx,2);
if i>l
    error('iԽ��');
end
if i~=1
    if bp(IntDNAx(i),IntDNAy(i))==0
        j=1;
        while j<=l-i && bp(IntDNAx(i+j),IntDNAy(i+j))
            IntValue=IntValue+1;
            j=j+1;
        end
    else
        return;
    end
else  %i==1
    if bp(IntDNAx(i),IntDNAy(i))==0
        j=1;
        while j<=l-i && bp(IntDNAx(i+j),IntDNAy(i+j))
            IntValue=IntValue+1;
            j=j+1;
        end
    else
        j=0;
        while j<=l-i && bp(IntDNAx(i+j),IntDNAy(i+j))
            IntValue=IntValue+1;
            j=j+1;
        end
    end
end
end
%------------------------------------------------------------------------------------------%


function Values= Continuity(DNAs)
% ��Shin2005,����DNAs��������,������ֵ������Values,������'-',tΪ�����Կ�����ֵ
global Con_T;
Con_T = 2;
[m,l]=size(DNAs);
Values=zeros(m,1);   % ��ʼ������ֵ������
t=Con_T;
for j=1:m
    Values(j)=0;
    for i=1:l-t+1
        temp=0;
        for a=0:3
            ca=c_a(DNAs(j,:),i,a);
            temp=temp+T(ca,t)^2;
        end
        Values(j)=Values(j)+temp;
    end
end
end

function IntValue = c_a( IntDNAx,i,a)
%Continuity�����и�������,�ӵ�i+1��λ�ÿ�ʼ�������ּ��a���Ӵ�����
%tΪ��ֵ����
IntValue=0;
l=size(IntDNAx,2);
if i>l
    error('iԽ��');
end
if i~=1
    if IntDNAx(i)~=a
        j=1;
        while j<=l-i && IntDNAx(i+j)==a
            IntValue=IntValue+1;
            j=j+1;
        end
    else
        return;
    end
else  %i==1
    if IntDNAx(i)~=a
        j=1;
        while j<=l-i && IntDNAx(i+j)==a
            IntValue=IntValue+1;
            j=j+1;
        end
    else
        j=0;
        while j<=l-i && IntDNAx(i+j)==a
            IntValue=IntValue+1;
            j=j+1;
        end
    end
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



%-----------------------------------�ظ��ԣ������ԣ�Լ��--------------------------------%
function con=continuity(x,numDNA)
c1=0;
c2=0;
cmax=0;
x(numDNA+1)=0;  %��֤ĩβ��ͬ
for i=1:numDNA
    if x(i)==x(i+1)
        c1=c1+1;
    else
        c2=c1;
        c1=0;
        if c2>=cmax
            cmax=c2;
        end
    end
end
con=cmax;
end
%----------------------------------------------------------------------------------%

%-----------------------------------GC����Լ��------------------------------------%
function gc=gccontent(x,numDNA)
p=0;
for i=1:20
    if x(i)=='G' || x(i)=='C'
        p=p+1;
    end
end
%p=(p/20)*100;
gc=p/20;
%gc=num2str(p)+"%";
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
                Values(k)=Values(k)+TT(Sigma_bp,pinlen(p,r,i,IntDNA)/ PINLEN);
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

function IntValue= TT(i,j)
%i>j����i�����򷵻�0
if i>j
    IntValue=i;
else
    IntValue=0;
end
end
%------------------------------------------------------------------------------------------%


%-----------------------����------------------------%
function betterdnaindex=evaluate(dna,numDNA,numSave)
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
%����һ�γ��ԡ�7�����������У�Tm40�֣�hanmingռ10�֡�simռ15�֡�H_mռ15�֡�conռ20�֣�����100��
%���ڶ��γ��ԡ�6�����������У�Tm40�֣�simռ20�֡�H_mռ20�֡�conռ20�֣�����100��
%�������γ��ԡ�6�����������У�Tm30�֣�simռ27�֡�H_mռ27�֡�conռ16�֣�����100��
%�����Ĵγ��ԡ�6�����������У�Tm38�֣�simռ21�֡�H_mռ21�֡�conռ20�֣�����100��
%������γ��ԡ�6�����������У�Tm20�֣�simռ30�֡�H_mռ30�֡�conռ20�֣�����100��
%�������γ��ԡ�6�����������У�Tm40�֣�simռ30�֡�H_mռ30�֡�conռ20�֣�����120��
%1��GC������50%��hairpin������0���ܷ�ֱ�Ӽ�Ϊ0��
%2��Tm����<=57��=40�֣���57~67��40���𼶵ݼ�����������ɢ�ݼ�������>=67��=0��
%4��sim������DNA��������0~5��=30����6~11��ÿλ��5,���ԣ���>=11���ܷ�ֱ�Ӽ�Ϊ0��
%5��H_m�����Գ���DNA��������0~5��=30����6~11��ÿλ��5,���ԣ���>=11���ܷ�ֱ�Ӽ�Ϊ0��
%6��con����=0��=20�֣���=1��=10�֣���=2��=5�֣���>=3��=0��
goal=[];
for i=1:numDNA
    goal(i)=0;
end
for i=1:numDNA
    %if hanming(dna(i,:),dna(:,:))<=20 && hanming(dna(i,:),dna(:,:))>15
    %   goal(i)=goal(i)+(hanming(dna(i,:))-15)*2;
    % else
    %    goal(i)=goal(i)+0;
    %end
    if continuity(dna(i,:),numDNA)==0
        goal(i)=goal(i)+20;
    elseif continuity(dna(i,:),numDNA)==1
        goal(i)=goal(i)+10;
    elseif continuity(dna(i,:),numDNA)==2
        goal(i)=goal(i)+5;
    else
        goal(i)=goal(i)+0;
    end
    if TmBioBox(dna(i,:))<=57
        goal(i)=goal(i)+40;
    elseif  TmBioBox(dna(i,:))>57&& TmBioBox(dna(i,:))<=67
        goal(i)=goal(i)+((67-TmBioBox(dna(i,:)))/10)*40;
    else
        goal(i)=goal(i)+0;
    end
    [sim,smax]=similarity(dna(i,:),dna(:,:),numDNA,i);
    [hme,hmax]=h_measure(dna(i,:),dna(:,:),numDNA);
    if sim/(numDNA-1)<=5
        goal(i)=goal(i)+30;
    elseif sim/(numDNA-1)>=6 && sim/(numDNA-1)<11
        goal(i)=goal(i)+(11-sim/(numDNA-1))*5;
    elseif sim/(numDNA-1)>=11 || smax>10
        goal(i)=0;
    end
    if hme/numDNA<=5
        goal(i)=goal(i)+30;
    elseif hme/numDNA>=6 && hme/numDNA<11
        goal(i)=goal(i)+(11-hme/numDNA)*5;
    elseif hme/numDNA>=11 || hmax>10 || sim/(numDNA-1)>=11 || smax>10
        goal(i)=0;
    end
    if gccontent(dna(i,:),numDNA)==0.5 && Hairpin(dna(i,:))==0
        goal(i)=goal(i);
    else
        goal(i)=0;
    end
    fprintf('��%d��DNA�ĵ÷�=%.2f\n',i,goal(i))
end
[sg,IX]=sort(goal,'descend');
[~,UIX]=unique(sg,'stable');
betterdnaindex=IX(UIX);
end
%----------------------------------------------------%


%����DNA�ܷ�
function GTotal=goaltotal(dna,numDNA)
goal=double(numDNA);
for i=1:numDNA
    goal(i)=0;
end
for i=1:numDNA
    %  if hanming(dna(i,:),dna(:,:))<=20 && hanming(dna(i,:),dna(:,:))>15
    % else
    %      goal(i)=goal(i)+0;
    %  end
    if continuity(dna(i,:),numDNA)==0
        goal(i)=goal(i)+20;
    elseif continuity(dna(i,:),numDNA)==1
        goal(i)=goal(i)+10;
    elseif continuity(dna(i,:),numDNA)==2
        goal(i)=goal(i)+5;
    else
        goal(i)=goal(i)+0;
    end
    if TmBioBox(dna(i,:))<=57
        goal(i)=goal(i)+40;
    elseif  TmBioBox(dna(i,:))>57&& TmBioBox(dna(i,:))<=67
        goal(i)=goal(i)+((67-TmBioBox(dna(i,:)))/10)*40;
    else
        goal(i)=goal(i)+0;
    end
    [sim,smax]=similarity(dna(i,:),dna(:,:),numDNA,i);
    [hme,hmax]=h_measure(dna(i,:),dna(:,:),numDNA);
    if sim/(numDNA-1)<=5
        goal(i)=goal(i)+30;
    elseif sim/(numDNA-1)>=6 && sim/(numDNA-1)<11
        goal(i)=goal(i)+(11-sim/(numDNA-1))*5;
    elseif sim/(numDNA-1)>=11 || smax>10
        goal(i)=0;
    end
    if hme/numDNA<=5
        goal(i)=goal(i)+30;
    elseif hme/numDNA>=6 && hme/numDNA<11
        goal(i)=goal(i)+(11-hme/numDNA)*5;
    elseif hme/numDNA>=11 || hmax>10 || sim/(numDNA-1)>=11 || smax>10
        goal(i)=0;
    end
    if gccontent(dna(i,:),numDNA)==0.5 && Hairpin(dna(i,:))==0
        goal(i)=goal(i);
    else
        goal(i)=0;
    end
end
GTotal=0;
for ii=1:numDNA
    GTotal=goal(ii)+GTotal;
end
end


%----------------------------------------����----------------------------------------------%
function d=Chemotaxis(dna,nc,numDNA)
dnaO=dna;
for i=1:nc   %��������
    fprintf('.')
    %-------------------�����������-------------------%
    for ii=1:numDNA
        % dnao=dna;
        % if hanming(dna(ii,:),dna(:,:))<20
        %     choose=ceil(rand(1,1)*20);
        %     if  dna(ii,choose)=='A'
        %         CH=['T','G','C'];
        %         change=CH(randperm(length(CH),1));
        %         dna(ii,choose)=change;
        %     elseif dna(ii,choose)=='T'
        %         CH=['A','G','C'];
        %         change=CH(randperm(length(CH),1));
        %         dna(ii,choose)=change;
        %     elseif dna(ii,choose)=='G'
        %         CH=['A','T','C'];
        %         change=CH(randperm(length(CH),1));
        %         dna(ii,choose)=change;
        %     elseif dna(ii,choose)=='C'
        %         CH=['A','T','G'];
        %        change=CH(randperm(length(CH),1));
        %         dna(ii,choose)=change;
        %    end
        % end
        % if hanming(dnao(ii,:),dnao(:,:))>hanming(dna(ii,:),dna(:,:))
        %     dna=dnao;   %����������º�������÷ָ�������
        % end
        %end
        %----------------------------------------------------%
        
        
        %--------------�����Ժ�H_measure����---------------%
        %for ii=1:20
        dnao=dna;
        [sim,smax]=similarity(dna(i,:),dna(:,:),numDNA,i);
        [hme,hmax]=h_measure(dna(i,:),dna(:,:),numDNA);
        [simo,smaxo]=similarity(dnao(i,:),dnao(:,:),numDNA,i);
        [hmeo,hmaxo]=h_measure(dnao(i,:),dnao(:,:),numDNA);
        if sim>0 || hme>0
            choose=ceil(rand(1,1)*20);%���ѡ��һλ
            if  dna(ii,choose)=='A'
                CH=['T','G','C'];
                change=CH(randperm(length(CH),1));
                dna(ii,choose)=change;
            elseif dna(ii,choose)=='T'
                CH=['A','G','C'];
                change=CH(randperm(length(CH),1));
                dna(ii,choose)=change;
            elseif dna(ii,choose)=='G'
                CH=['A','T','C'];
                change=CH(randperm(length(CH),1));
                dna(ii,choose)=change;
            elseif dna(ii,choose)=='C'
                CH=['A','T','G'];
                change=CH(randperm(length(CH),1));
                dna(ii,choose)=change;
            end
        end
        if sim>simo || hme>hmeo
            dna=dnao;   %����������������Ե÷ָ�������   �����������H_measure�÷ָ�������
        elseif smax>smaxo || hmax>hmaxo
            dna=dnao;
        elseif smax>=10 || hmax>=10
            dna=dnao;
        end
        % end
        %-------------------------------------------------%
        
        %----------------�ظ���(������)����----------------%
        %  for ii=1:20
        %dnao=dna;
        if continuity(dna(ii,:),numDNA)>0
            for iii=1:19
                if dna(ii,iii)==dna(ii,iii+1) && dna(ii,iii)=='A'
                    CH=['T','G','C'];
                    change=CH(randperm(length(CH),1));
                    dna(ii,iii+1)=change;
                elseif dna(ii,iii)==dna(ii,iii+1) && dna(ii,iii)=='T'
                    CH=['A','G','C'];
                    change=CH(randperm(length(CH),1));
                    dna(ii,iii+1)=change;
                elseif dna(ii,iii)==dna(ii,iii+1) && dna(ii,iii)=='G'
                    CH=['A','T','C'];
                    change=CH(randperm(length(CH),1));
                    dna(ii,iii+1)=change;
                elseif dna(ii,iii)==dna(ii,iii+1) && dna(ii,iii)=='C'
                    CH=['A','T','G'];
                    change=CH(randperm(length(CH),1));
                    dna(ii,iii+1)=change;
                end
                if continuity(dna(ii,:),numDNA)>continuity(dnao(ii,:),numDNA)
                    dna=dnao;   %����������������Ե÷ָ�������
                end
            end
        end
        % end
        %-------------------------------------------------%
        
        %---------------GC��������-------------%
        %for ii=1:20
        if gccontent(dna(ii,:),numDNA)<0.5
            r=0;
            remark=[];
            for iii=1:20
                if dna(ii,iii)=='A' || dna(ii,iii)=='T'
                    r=r+1;
                    remark(r)=iii;
                end
            end
            n=length(remark);
            choose=ceil(rand(1,1)*n);%���ѡ��һλ
            CH=['C','G'];
            change=CH(randperm(length(CH),1));%����޸ļ��
            dna(ii,choose)=change;
        end
        if gccontent(dna(ii,:),numDNA)>0.5
            r=0;
            remark=[];
            for iii=1:20
                if dna(ii,iii)=='G' || dna(ii,iii)=='C'
                    r=r+1;
                    remark(r)=iii;
                end
            end
            n=length(remark);
            choose=ceil(rand(1,1)*n);%���ѡ��һλ
            CH=['A','T'];
            change=CH(randperm(length(CH),1));%����޸ļ��
            dna(ii,choose)=change;
        end
        if gccontent(dna(ii,:),numDNA)==0.5
        end
        %----------------------------------------------%
        
        GOT=goaltotal(dnaO,numDNA);
        GT=goaltotal(dna,numDNA);
        
        if GOT>GT
            dna=dnaO;
        else
            dnaO=dna;
        end
    end
end
fprintf('\n')

fprintf('\n')
d=dna;
end

%----------------------------------------------����--------------------------------------------%
function [re,ped]=reproduce(dna)
re=[dna;dna];
ped=rand(1,1);
end
%-----------------------------------------------------------------------------------------------%

%-------------------------------------------�������--------------------------------------------%
function e=eliminate(dna)
[m,~]=size(dna);
i=randi([m/4,3*m/4],1,m/4);            %���ѡDNA����m�����ķ�֮һm���ķ�֮��m֮����ķ�֮һm��
choose=ceil(rand(1,2)*20);%���ѡ���λ
CH=['A','T','G','C'];
change=CH(randperm(length(CH),1));
dna(i,choose)=change;
e=dna;
end
%------------------------------------------------------------------------------------------------%