function IntValue = ceq( IntDNAx,IntDNAy,i)
%�ӵ�i��λ�ÿ�ʼ����Ϊc��IntDNAx��IntDNAy�Ĺ����Ӵ��ĳ���;IN:IntDNAx��IntDNAy,��������DNA,i,��ʼλ��;OUT:IntValue,������ִ�����
IntValue=0;
l=size(IntDNAx,2);
if i>l
    error('iԽ��');
end
if i~=1
    if eqq(IntDNAx(i),IntDNAy(i))==0
        j=1;
        while j<=l-i && eqq(IntDNAx(i+j),IntDNAy(i+j))
            IntValue=IntValue+1;
            j=j+1;
        end
    else
        return;
    end
else %i==1
    if eqq(IntDNAx(i),IntDNAy(i))==0
        j=1;
        while j<=l-i && eqq(IntDNAx(i+j),IntDNAy(i+j))
            IntValue=IntValue+1;
            j=j+1;
        end        
    else
        j=0;
        while j<=l-i && eqq(IntDNAx(i+j),IntDNAy(i+j))
            IntValue=IntValue+1;
            j=j+1;
        end
    end    
end
end


