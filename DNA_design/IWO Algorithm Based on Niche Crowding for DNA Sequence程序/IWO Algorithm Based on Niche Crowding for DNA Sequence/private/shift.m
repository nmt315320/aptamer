function IntDNAShifted= shift( IntDNA_nb,i )
%将IntDNA_nb平移i位，i>0右移，i<0,左移;空位用'-'(5)添加
if i==0
    IntDNAShifted=IntDNA_nb;
    return;
end
l=size(IntDNA_nb,2);%DNA长度
temp=zeros(1,abs(i))+5;
if (0<i) && (i<l)   %右移i位添加i个'-'
    IntDNAShifted=[temp, IntDNA_nb(1:l-i)];
    return;
end
if (i<0) && (i>-l)  %左移i位添加i个'-'
    IntDNAShifted=[IntDNA_nb(abs(i)+1:l),temp];
    return;
end
if abs(i)>=l %全部替换为'-'
    IntDNAShifted=zeros(1,l)+5;
    return;
end
end

