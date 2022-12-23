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

