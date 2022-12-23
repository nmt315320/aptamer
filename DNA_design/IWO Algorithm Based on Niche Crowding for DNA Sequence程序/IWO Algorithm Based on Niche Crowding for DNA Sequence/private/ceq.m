function IntValue = ceq( IntDNAx,IntDNAy,i)
%从第i个位置开始长度为c的IntDNAx与IntDNAy的公共子串的长度;IN:IntDNAx、IntDNAy,整数编码DNA,i,起始位置;OUT:IntValue,最长公共字串长度
IntValue=0;
l=size(IntDNAx,2);
if i>l
    error('i越界');
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


