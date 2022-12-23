function Values= Continuity(DNAs)
% 依Shin2005,计算DNAs的连续性,返回数值列向量Values,不包括'-',t为连续性控制阈值
global Con_T;
Con_T = 2;
[m,l]=size(DNAs);
Values=zeros(m,1);   % 初始化连续值列向量
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
%Continuity函数中辅助函数,从第i+1个位置开始连续出现碱基a的子串长度
%t为阈值控制
IntValue=0;
l=size(IntDNAx,2);
if i>l
    error('i越界');
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
%i>j返回i，否则返回0
if i>j
    IntValue=i;
else
    IntValue=0;
end
end

