%ʹ��˵����������4.���ղ�ͬ��numSave����������ճ��                                                           %
%                         5.������У���Ȼ�ȴ������д�����ʾ������ʱ��                                             %
%                            �мǣ�����������ֻ��������ʾ��Ž�չ�ģ����Ը��ݽ�������Ϣȥ�ȵ�裻    %
%                            �мǣ�������������ʱ��󣬹رս��������ڣ�                                             %
%                         6.��¼���ÿһ����  ����š�  ��hme��  ��sim��  �����ĺϼ���hme��sim                %
%                            hme��sim�������е����ݱ�ʾ����DNA���бȽϹ����е����ֵ                    %  
%                         7.�õ�ǰ����е���ţ�ȥBFOA��ǰһ���ӳ��򣩵Ľ����Ѱ����ͬ��ź�����   %
%                            �е���š���ϸ����仰��                                                                          %
%                         8.������Щ�����е���ţ�ת�����һ�θ��ƵĽ�����棬ȥ���Ӧ�ĺ��ĸ�Լ�� %
%                            Continuity��GC_content��Tm��Hairpinnum
%                            �мǣ������������һ�θ��ƵĽ����ֻ�����ĸ�Լ��������                          %
%                                                       ������ˣ������򵥲�                                                      %





tic
DNA1_2=['CACGACTCTAGCGATCATAG';'CTGCGAGAGTAGATAGACAG'];
DNA3_4=['CATAGCTAGCAGTACTGTGC';'CGTCTATCTGCTAGCATGAC'];
DNA5_6=['GTCATACGAGACTCTGAGTG';'GACACGTCAGATACATCTGC'];
DNA7_8=['GTCATACGTGAGTCTGAGTG';'CGCTATCTGATAGCGACTAG'];
DNA9_10=['TAGCAGCGATGCTACTAGAC';'GATACGTGAGACACATCTGC'];
DNA11_12=['';''];
DNA13_14=['';''];
DNA15_16=['';''];
DNA17_18=['';''];
DNA19_20=['';''];
DNA21_22=['';''];
DNA23_24=['';''];
DNA25_26=['';''];
DNA27_28=['';''];
DNA29_30=['';''];
DNA31_32=['';''];
DNA33_34=['';''];
DNA35_36=['';''];
DNA37_38=['';''];
DNA39_40=['';''];
DNA1=[DNA1_2;DNA3_4;DNA5_6;DNA7_8;DNA9_10];
DNA2=[DNA11_12;DNA13_14;DNA15_16;DNA17_18;DNA19_20];
DNA3=[DNA21_22;DNA23_24;DNA25_26;DNA27_28;DNA29_30];
DNA4=[DNA31_32;DNA33_34;DNA35_36;DNA37_38;DNA39_40];
DNA=[DNA1;DNA2;DNA3;DNA4];
[DNAnum,~]=size(DNA);


d1='ACACCAGCACACCAGAAACA';
d2='GTTCAATCGCCTCTCGGTAT';
d3='GCTACCTCTTCCACCATTCT';
d4='GAATCAATGGCGGTCAGAAG';
d5='TTGGTCCGGTTATTCCTTCG';
d6='CCATCTTCCGTACTTCACTG';
d7='TTCGACTCGGTTCCTTGCTA';
dnayang=[d1;d2;d3;d4;d5;d6;d7];

ds1='CTCTTCATCCACCTCTTCTC';
ds2='CTCTCATCTCTCCGTTCTTC';
ds3='TATCCTGTGGTGTCCTTCCT';
ds4='ATTCTGTTCCGTTGCGTGTC';
ds5='TCTCTTACGTTGGTTGGCTG';
ds6='GTATTCCAAGCGTCCGTGTT';
ds7='AAACCTCCACCAACACACCA';
dnashin=[ds1;ds2;ds3;ds4;ds5;ds6;ds7];

dy1='GATGGATTTACCTTGCACCT';
dy2='CCTTCTCTCGTCTTCATACA';
dy3='ACGATCGATTAATGGGAGTC';
dy4='ATAAGTAGGGACTGCTCTAC';
dy5='CCTAAGAACACAGGGCATAG';
dy6='CTGGAAGCGTTTGCTAACTT';
dy7='GCAGATTCCCGGATACTCAG';
dnay=[dy1;dy2;dy3;dy4;dy5;dy6;dy7];


dn=7;
finalchooseDNAanalysis(DNA,DNAnum);
%����ע�͵����ǶԱ�ʵ����%
%finalchooseDNAanalysis(dnayang,dn);
%finalchooseDNAanalysis(dnashin,dn);
%finalchooseDNAanalysis(dnay,dn);

toc

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

function finalchooseDNAanalysis(dna,dnanum)
so=0;
ho=0;
record=[];
bestdna=[];
choose=nchoosek((1:dnanum),7);
w=waitbar(0,'��ʼ'); 
display=factorial(dnanum)/(factorial(dnanum-7)*factorial(7));
for i=1:factorial(dnanum)/(factorial(dnanum-7)*factorial(7))
    if i/display<0.9
    waitbar(i/display,w,['ū���Ѿ�����' num2str(i) '���ȶ��ˣ��������ٺȵ��']);
    elseif i/display>=0.9
    waitbar(i/display,w,['ū���Ѿ�����' num2str(i) '���ȶ��ˣ���Ҫ���ˣ�']);   
    elseif i==diaplay
    waitbar(1,w,'�������ӣ����ˣ���');
    end
    Choose=choose(i,:);
    choosedna=dna(Choose,:);
    choosednanum=7;
    s=0;
    h=0;
    for choosei=1:choosednanum
        [sim,~]=similarity(choosedna(choosei,:),choosedna(:,:),choosednanum,choosei);
        s=s+sim;
        [hme,~]=h_measure(choosedna(choosei,:),choosedna(:,:),choosednanum);
        h=h+hme;
    end
    if i==1
        so=s;
        ho=h;
    end
    if s<=so && h<=ho
        so=s;
        ho=h;
        record=Choose;
        bestdna=dna(record,:);
    end
end
stotal=0;
htotal=0;
for chooseii=1:7
    [Ssim,Smax]=similarity(bestdna(chooseii,:),bestdna(:,:),choosednanum,chooseii);
    stotal=Ssim+stotal;
    [Hhme,Hmax]=h_measure(bestdna(chooseii,:),bestdna(:,:),choosednanum);
    htotal=Hhme+htotal;
    fprintf('  %s  ��%d����(���)hme��%d(%d)    sim��%d(%d)\n',bestdna(chooseii,:),record(chooseii),Hhme,Hmax,Ssim,Smax)
end
fprintf('�ϼ��ܵ�hme��%d    sim��%d\n',htotal,stotal)
end