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

