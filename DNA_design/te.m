clc
d1='CTAGAGAGAGAGAGAGAGAG';
d2='CTAGAGAGAGAGAGACAGAG';
num=2;
similarity(d1,d2,num)


function S=similarity(xo,yo,numDNA)
s=[];
for i=1:40
    s(i)=0;
end

for k=1:19  %k表示y向'右'移动距离
    x=xo;
    y=circshift(yo,k);
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
        if  x(i)==y(i)
            s(19+kk)=s(kk+19)+1;
        end
    end
end
for i=1:20
    if  xo(i)==yo(i)
        s(40)=s(40)+1;
    end
end
s
[S,~]=max(s);
end