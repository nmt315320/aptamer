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

