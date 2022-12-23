clc
d1='CTAGAGAGAGAGAGAGAGAG';
d2='CTAGAGAGAGAGAGACAGAG';
num=2;
H_measure(d1,[d1;d2],num)


function hm=H_measure(xo,yo,numDNA)
s=[];
for i=1:40
    s(i)=0;
end
htotal=0;
for ci=1:numDNA   %对y取补
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
yo
for ii=1:numDNA
    for k=1:19  %k表示y向'右'移动距离
        x=xo;
        y=circshift(yo(ii,:),k);
        for i=k:20
            if  x(i)==y(i)
                s(k)=s(k)+1;
            end
        end
    end
    for kk=1:19 %kk表示y向'左'移动距离(用x向右移动替代y向左移动)
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
    [ss,~]=max(s);%hmp表示出现最大相似的位置（20之前为y右移，21之后为y左移）
    s(:,:)=0;
    htotal=htotal+ss;
end
hm=htotal;
end